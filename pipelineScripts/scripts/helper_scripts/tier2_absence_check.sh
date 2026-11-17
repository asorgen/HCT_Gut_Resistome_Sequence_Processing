#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=tier2_absence_check

# Sequence-based whole-sample absence test for tier-2 ARG gain candidates
# (bug-to-bug ARG transfer, module 4.15 in HCT_Gut_Resistome_Species_Attribution).
#
# 4.15 calls a gene absent from a sample by comparing RGI names. RGI reports only the closest
# CARD allele, so a sister allele of the same gene reads as absence and the "gain" is an
# artefact: CfxA2 was scored as gained in D21276 while CfxA3 sat in the earlier sample at 100%
# identity. Grouping beta-lactamases by name stem catches some of it, but not tet(X) vs tet(X1)
# or the tet(O/32/O) mosaics, since elsewhere a shared stem names genuinely different genes
# (sul1 / sul2, dfrA1 / dfrA12, tet(40) / tet(41)).
#
# This searches for the gene's own sequence instead. For each candidate (a gene-free sample in a
# tier-2 trajectory) the query is the gene's ORF on the carrying contig named in the candidate
# file, and the target is the candidate sample's whole assembly (1.2), so the result does not
# depend on binning or on what RGI chose to call anything.
#
# Queries are batched by target sample: one minimap2 run per assembly, not per candidate, since
# loading the assembly dominates the cost.
#
# Input : <4.15 output>/Tier2_absence_candidates.tsv, written by 4.15. Columns:
#         Cohort, PatientID, gtdb_species, ARO, Gene, SampleID (the sample to search),
#         Day, SampleID_query, contig_query (where to take the gene sequence from).
# Output: <4.15 output>/Tier2_gene_absence.tsv, one row per candidate:
#         the key columns plus gene_cov_before (share of the gene's length aligned),
#         gene_id_before (identity of the best hit, 0 when there is none),
#         gene_hit_contig and gene_hit_rgi (the RGI calls on the hit contig, which name the
#         sister allele responsible), and query_len.
#         4.15 reads it back and treats the gene as present at cov >= min_gene_cov (0.5) and
#         identity >= min_gene_id (0.95).
#
# Usage (submits itself; run from the login node):
#   bash tier2_absence_check.sh <log_file> <path/to/Tier2_absence_candidates.tsv>

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: tier2_absence_check.sh <log_file> <Tier2_absence_candidates.tsv>}
    cand=$(realpath "${2:?Usage: tier2_absence_check.sh <log_file> <Tier2_absence_candidates.tsv>}")
    sbatch -o "$log" --export=ALL,AC_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",AC_CAND="${cand}" \
        "${SCRIPT_DIR}/tier2_absence_check.sh"
    exit $?
fi
CONFIG_DIR=${AC_CONFIG_DIR:?}
CAND=${AC_CAND:?}
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
module load minimap2/2.28
echo "minimap2 $(minimap2 --version)"

cohorts=()
for c in Duke_short UNC_short Heston_short; do
    v=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$datasetDir|$evaluationDir")
    cohorts+=("${c%_short}=$v")
done

work=${TMPDIR:-/tmp}/ac_${SLURM_JOB_ID}
mkdir -p "$work"; trap 'rm -rf "$work"' EXIT

$PY - "$CAND" "$work" "$THREADS" "${cohorts[@]}" <<'PYEOF'
import csv, os, re, subprocess, sys
from collections import defaultdict
cand_tsv, work, threads = sys.argv[1:4]
coh = {}
for a in sys.argv[4:]:
    k, v = a.split("=", 1)
    ds, ev = v.split("|")
    coh[k] = {"ds": ds, "eval": ev}
out = os.path.join(os.path.dirname(cand_tsv), "Tier2_gene_absence.tsv")
MIN_ALN = 50          # ignore scraps; the genes here are 500-2000 bp

def cohort_of(row, sid):
    return "Heston" if row["Cohort"] == "Heston" else ("UNC" if sid.startswith("BMT") else "Duke")

def asm_fa(c, sid):
    return os.path.join(coh[c]["eval"], f"{sid}_final_assembly.fasta")

def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

rows = read_tsv(cand_tsv)
print(f"{len(rows)} candidates, {len({r['SampleID'] for r in rows})} distinct samples to search")

# ---- query sequences: the gene's ORF on the named carrying contig ----
# Prodigal ORF coordinates come from 5.5 <ID>_genes.faa; the gene is the RGI hit(s) for that ARO
# on that contig at >= 80% identity, the same floor 4.15 uses.
need = defaultdict(set)                       # (cohort, query sample) -> {(contig, ARO, Gene)}
for r in rows:
    cq = cohort_of(r, r["SampleID_query"])
    need[(cq, r["SampleID_query"])].add((r["contig_query"], r["ARO"], r["Gene"]))

qseq = {}                                     # (ARO, Gene, query sample, contig) -> sequence
for (c, sid), want in sorted(need.items()):
    contigs = {w[0] for w in want}
    seqs, cur = {}, None
    with open(asm_fa(c, sid)) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                cur = name if name in contigs else None
                if cur:
                    seqs[cur] = []
            elif cur:
                seqs[cur].append(line.strip())
    seqs = {k: "".join(v) for k, v in seqs.items()}
    d55 = os.path.join(coh[c]["ds"], "5.5_AA_amr_assembly", sid)
    orfs = {}
    with open(os.path.join(d55, f"{sid}_genes.faa")) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].split(" # ")
                oid = f[0]
                if oid.rsplit("_", 1)[0] in contigs:
                    orfs[oid] = (int(f[1]), int(f[2]))
    spans = defaultdict(list)                 # (contig, Best_Hit_ARO) -> [(start, end)]
    for r in read_tsv(os.path.join(d55, f"{sid}.rgi.txt")):
        oid = r["ORF_ID"].split(" # ")[0]
        if oid in orfs and float(r["Best_Identities"] or 0) >= 80:
            spans[(oid.rsplit("_", 1)[0], r["Best_Hit_ARO"])].append(orfs[oid])
    for contig, aro, gene in sorted(want):
        sp = spans.get((contig, gene))
        if not sp or contig not in seqs:
            print(f"  WARN no >=80% {gene} ORF on {sid} {contig}; candidates using it get no result")
            continue
        s, e = min(x[0] for x in sp) - 1, max(x[1] for x in sp)
        qseq[(aro, gene, sid, contig)] = seqs[contig][s:e]

# ---- one minimap2 per target sample, all its queries at once ----
by_target = defaultdict(list)
for i, r in enumerate(rows):
    by_target[(cohort_of(r, r["SampleID"]), r["SampleID"])].append(i)

res = {}
rgi_cache = {}
def rgi_of(c, sid):
    if (c, sid) not in rgi_cache:
        d = defaultdict(list)
        for r in read_tsv(os.path.join(coh[c]["ds"], "5.5_AA_amr_assembly", sid, f"{sid}.rgi.txt")):
            d[r["ORF_ID"].split(" # ")[0].rsplit("_", 1)[0]].append(f"{r['Best_Hit_ARO']}({r['Best_Identities']})")
        rgi_cache[(c, sid)] = d
    return rgi_cache[(c, sid)]

def covered(iv, n):
    tot, cs, ce = 0, None, None
    for s, e in sorted(iv):
        if ce is None or s > ce:
            if ce is not None:
                tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    return (tot + (ce - cs)) if ce is not None else 0

for n, ((c, sid), idx) in enumerate(sorted(by_target.items()), 1):
    qs = {}
    for i in idx:
        r = rows[i]
        k = (r["ARO"], r["Gene"], r["SampleID_query"], r["contig_query"])
        if k in qseq:
            qs[f"{r['ARO']}|{r['Gene']}"] = qseq[k]
    if not qs:
        continue
    qf = os.path.join(work, "q.fa")
    with open(qf, "w") as fo:
        for name, s in qs.items():
            fo.write(f">{name}\n{s}\n")
    p = subprocess.run(["minimap2", "-x", "asm20", "-c", "-t", threads, asm_fa(c, sid), qf],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"minimap2 failed on {sid}: {p.stderr[-300:]}")
    hits = defaultdict(list)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if int(f[10]) >= MIN_ALN:
            hits[f[0]].append(f)
    for i in idx:
        r = rows[i]
        key = f"{r['ARO']}|{r['Gene']}"
        if key not in qs:
            continue
        h = hits.get(key, [])
        qlen = len(qs[key])
        if not h:
            res[i] = ("0.00", "0.000", "", "", qlen)
        else:
            best = max(h, key=lambda x: int(x[10]))
            tg = sorted({x[5] for x in h})
            rg = rgi_of(c, sid)
            res[i] = (f"{covered([(int(x[2]), int(x[3])) for x in h], qlen) / qlen:.2f}",
                      f"{int(best[9]) / int(best[10]):.3f}",
                      ",".join(tg)[:300],
                      "; ".join(", ".join(rg.get(t, ["none"])) for t in tg)[:300],
                      qlen)
    if n % 25 == 0 or n == len(by_target):
        print(f"  {n}/{len(by_target)} samples searched", flush=True)

cols = list(rows[0]) + ["gene_cov_before", "gene_id_before", "gene_hit_contig", "gene_hit_rgi", "query_len"]
with open(out, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(cols)
    n_written = 0
    for i, r in enumerate(rows):
        if i not in res:
            continue
        w.writerow([r[c] for c in rows[0]] + list(res[i]))
        n_written += 1

present = sum(1 for v in res.values() if float(v[0]) >= 0.5 and float(v[1]) >= 0.95)
print(f"\n{out}: {n_written} of {len(rows)} candidates checked")
print(f"  gene already present (cov >= 0.5, id >= 0.95): {present}  <- these are NOT absences")
print(f"  gene absent: {n_written - present}")
if n_written < len(rows):
    print(f"  {len(rows) - n_written} candidates had no usable query ORF (see WARNs above); 4.15 keeps its own test for them")
# Strip the trailing "(identity)" only - gene names contain parentheses, e.g. tet(X1), aad(6).
nm = sorted({re.sub(r"\([0-9.]+\)$", "", x) for v in res.values()
             if float(v[0]) >= 0.5 and float(v[1]) >= 0.95 and v[3]
             for x in v[3].replace("; ", ", ").split(", ")})
print(f"  RGI names found at the present sites (sister alleles): {', '.join(nm[:25])}")
PYEOF
echo "done"
