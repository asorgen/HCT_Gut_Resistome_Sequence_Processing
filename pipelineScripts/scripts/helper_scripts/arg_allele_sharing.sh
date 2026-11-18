#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=arg_allele_sharing

# Cross-family allele sharing of ARGs, by comparing the gene copies themselves
# (bug-to-bug ARG transfer, module 4.15 in HCT_Gut_Resistome_Species_Attribution).
#
# 4.15 builds its shortlist from the 4.14 candidates whose allele is shared across >= 3 bacterial
# families (each carrying it in >= 2 patients). Locally it can only approximate "same allele" by
# the identity to the CARD reference: two copies 98.5% from the reference may still differ from
# each other. This replaces the proxy with the direct test, following Smillie et al. 2011 (Nature
# 480:241; recent transfer = >= 99% identity over >= 500 bp between distantly related genomes):
#   1. Pull every listed copy's ORF (nucleotide, from the assembly, oriented by Prodigal's strand).
#   2. Per gene, collapse identical sequences, then align all unique sequences against each other
#      (minimap2 -x asm20, all-vs-all).
#   3. Link two copies at >= MIN_ID identity over >= min(MIN_LEN, MIN_FRAC x the shorter copy)
#      (the length floor is relaxed for genes shorter than 500 bp, e.g. dfrF ~500 bp).
#   4. Alleles = connected components of those links (single linkage).
#   5. Per allele, count GTDB families (suffix stripped: Lachnospiraceae_A -> Lachnospiraceae) that
#      carry it in >= MIN_PAT patients; families_same_allele = the most for any allele of the gene.
#
# Input : <4.15 output>/Allele_sharing_queries.tsv, written by 4.15: one row per ARG hit on a good
#         MAG (>= 70% complete, <= 5% contaminated) for the 4.14 gene set. Columns: Source
#         (Duke / UNC / Heston), Cohort, PatientID, SampleID, MAG, contig, ORF_ID (Prodigal ID on
#         the contig), ARO, Gene, Best_Identities, gtdb_family, gtdb_genus, gtdb_species.
# Sources: <evaluationDir>/<SampleID>_final_assembly.fasta (contigs) and
#         <datasetDir>/5.5_AA_amr_assembly/<SampleID>/<SampleID>_genes.faa (ORF coordinates and strand),
#         as in tier2_absence_check.sh.
# Output: next to the input
#   Allele_sharing.tsv         per gene: ARO, Gene, n_copies, n_unique_seqs, n_alleles,
#                              families_same_allele, genera_same_allele, top_allele_families,
#                              n_linked_pairs. 4.15 reads it back in place of the proxy.
#   Allele_sharing_copies.tsv  per copy: the query columns plus allele (component id), seq_len.
#
# Usage (submits itself; run from the login node):
#   bash arg_allele_sharing.sh <log_file> <path/to/Allele_sharing_queries.tsv>

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: arg_allele_sharing.sh <log_file> <Allele_sharing_queries.tsv>}
    q=$(realpath "${2:?Usage: arg_allele_sharing.sh <log_file> <Allele_sharing_queries.tsv>}")
    sbatch -o "$log" --export=ALL,AS_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",AS_QUERIES="${q}" \
        "${SCRIPT_DIR}/arg_allele_sharing.sh"
    exit $?
fi
CONFIG_DIR=${AS_CONFIG_DIR:?}
QUERIES=${AS_QUERIES:?}
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
module load minimap2/2.28
echo "minimap2 $(minimap2 --version)"
echo "queries: $QUERIES"

cohorts=()
for c in Duke_short UNC_short Heston_short; do
    v=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$datasetDir|$evaluationDir")
    cohorts+=("${c%_short}=$v")
done

work=${TMPDIR:-/tmp}/as_${SLURM_JOB_ID}
mkdir -p "$work"; trap 'rm -rf "$work"' EXIT

$PY - "$QUERIES" "$work" "$THREADS" "${cohorts[@]}" <<'PYEOF'
import csv, os, re, subprocess, sys
from collections import defaultdict
q_tsv, work, threads = sys.argv[1:4]
coh = {}
for a in sys.argv[4:]:
    k, v = a.split("=", 1)
    ds, ev = v.split("|")
    coh[k] = {"ds": ds, "eval": ev}
MIN_ID, MIN_LEN, MIN_FRAC, MIN_PAT = 0.99, 500, 0.9, 2
outdir = os.path.dirname(q_tsv)

def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s):
    return s.translate(COMP)[::-1]

def strip_sfx(x):
    return re.sub(r"_[A-Z]+$", "", x or "")

rows = read_tsv(q_tsv)
cols = list(rows[0].keys())
print(f"{len(rows)} gene copies, {len({r['ARO'] for r in rows})} genes, {len({r['SampleID'] for r in rows})} samples")

# ---- 1. extract each copy's ORF ----
by_sample = defaultdict(list)
for i, r in enumerate(rows):
    by_sample[(r["Source"], r["SampleID"])].append(i)
seq = {}
missing = 0
for n, ((src, sid), idx) in enumerate(sorted(by_sample.items()), 1):
    want_orf = {rows[i]["ORF_ID"] for i in idx}
    want_ctg = {rows[i]["contig"] for i in idx}
    orfs = {}
    faa = os.path.join(coh[src]["ds"], "5.5_AA_amr_assembly", sid, f"{sid}_genes.faa")
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].split(" # ")
                if f[0] in want_orf:
                    orfs[f[0]] = (int(f[1]), int(f[2]), int(f[3]))
    contigs, cur = {}, None
    with open(os.path.join(coh[src]["eval"], f"{sid}_final_assembly.fasta")) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                cur = name if name in want_ctg else None
                if cur:
                    contigs[cur] = []
            elif cur:
                contigs[cur].append(line.strip())
    contigs = {k: "".join(v) for k, v in contigs.items()}
    for i in idx:
        r = rows[i]
        o, c = orfs.get(r["ORF_ID"]), contigs.get(r["contig"])
        if not o or c is None:
            missing += 1
            continue
        s = c[o[0] - 1:o[1]]
        seq[i] = (revcomp(s) if o[2] == -1 else s).upper()
    if n % 100 == 0:
        print(f"  extracted {n}/{len(by_sample)} samples", flush=True)
print(f"ORFs extracted: {len(seq)} of {len(rows)} ({missing} not found)")

# ---- 2-4. per gene: dedupe, all-vs-all, link, components ----
class UF:
    def __init__(s, n): s.p = list(range(n))
    def f(s, x):
        while s.p[x] != x:
            s.p[x] = s.p[s.p[x]]; x = s.p[x]
        return x
    def u(s, a, b): s.p[s.f(a)] = s.f(b)

by_gene = defaultdict(list)
for i in seq:
    by_gene[rows[i]["ARO"]].append(i)
allele_of, summary = {}, []
for aro, idx in sorted(by_gene.items()):
    uniq = {}
    for i in idx:
        uniq.setdefault(seq[i], len(uniq))
    useqs = sorted(uniq, key=uniq.get)
    uf = UF(len(useqs))
    n_links = 0
    if len(useqs) > 1:
        fa = os.path.join(work, "g.fa")
        with open(fa, "w") as fo:
            for k, s in enumerate(useqs):
                fo.write(f">{k}\n{s}\n")
        p = subprocess.run(["minimap2", "-x", "asm20", "-c", "-X", "-N", "100000", "-p", "0", "-t", threads, fa, fa],
                           capture_output=True, text=True)
        if p.returncode != 0:
            sys.exit(f"minimap2 failed on {aro}: {p.stderr[-300:]}")
        for line in p.stdout.splitlines():
            f = line.split("\t")
            qn, ql, tn, tl, match, blk = int(f[0]), int(f[1]), int(f[5]), int(f[6]), int(f[9]), int(f[10])
            if blk and match / blk >= MIN_ID and blk >= min(MIN_LEN, MIN_FRAC * min(ql, tl)):
                uf.u(qn, tn); n_links += 1
    comp = {}
    for i in idx:
        allele_of[i] = f"{aro}_a{uf.f(uniq[seq[i]])}"
        comp.setdefault(allele_of[i], []).append(i)
    best = (0, 0, "")
    for a, members in comp.items():
        fam_pat, gen_pat = defaultdict(set), defaultdict(set)
        for i in members:
            r = rows[i]
            pt = (r["Cohort"], r["PatientID"])
            fam_pat[strip_sfx(r["gtdb_family"])].add(pt)
            gen_pat[strip_sfx((r["gtdb_genus"] or "").strip())].add(pt)
        fams = sorted(f for f, ps in fam_pat.items() if f and len(ps) >= MIN_PAT)
        gens = [g for g, ps in gen_pat.items() if g and len(ps) >= MIN_PAT]
        if (len(fams), len(gens)) > best[:2]:
            best = (len(fams), len(gens), "; ".join(fams))
    summary.append([aro, rows[idx[0]]["Gene"], len(idx), len(useqs), len(comp), best[0], best[1], best[2], n_links])
    print(f"  {rows[idx[0]]['Gene'][:40]:<40} copies {len(idx):>5}  unique {len(useqs):>4}  alleles {len(comp):>4}  "
          f"families_same_allele {best[0]:>2}  genera {best[1]:>2}", flush=True)

with open(os.path.join(outdir, "Allele_sharing.tsv"), "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(["ARO", "Gene", "n_copies", "n_unique_seqs", "n_alleles", "families_same_allele", "genera_same_allele",
                "top_allele_families", "n_linked_pairs"])
    w.writerows(sorted(summary, key=lambda x: -x[5]))
with open(os.path.join(outdir, "Allele_sharing_copies.tsv"), "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(cols + ["allele", "seq_len"])
    for i, r in enumerate(rows):
        w.writerow([r[c] for c in cols] + [allele_of.get(i, "NA"), len(seq[i]) if i in seq else "NA"])
print(f"\nwrote {outdir}/Allele_sharing.tsv ({len(summary)} genes) and Allele_sharing_copies.tsv")
PYEOF
echo "done"
