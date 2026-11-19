#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=donor_contig_aln

# Donor-vs-recipient contig alignment for the 4.17 donor-present gains
# (bug-to-bug ARG transfer, HCT_Gut_Resistome_Species_Attribution).
#
# The ANI strain check cannot tell a real transfer between two resident genomes from one mobile
# element binned to different species on different days: BMT104 vanB passes both ANI tools and is
# a known rebinning artefact. The ANI test asks whether the RECIPIENT is the same strain over
# time; this asks whether the DONOR and RECIPIENT copies are the same physical sequence.
#
# For each event x donor copy, the donor's contig is aligned to the recipient's contig
# (minimap2 -x asm5) and the result is described by:
#   * aligned fraction of each contig, and of the SHORTER one (the discriminating number when the
#     two contigs differ greatly in length, e.g. H37: 1.9 kb donor vs 10.8 kb recipient - there
#     "one element binned twice" shows as the shorter contig aligning almost entirely, not as a
#     near-full match of both);
#   * the longest single aligned block, and whether each copy's ARG ORF falls inside it;
#   * how much sequence lies OUTSIDE that block on each contig (the flanks). Divergent flanks on
#     both sides are what distinguishes two genomes sharing an element from one sequence counted
#     twice.
#
# Verdicts (thresholds are the named parameters below). The first three all mean "the shorter
# contig is >= MIN_CONTAINED aligned into the longer one", i.e. the two samples hold ONE sequence
# rather than two genomes' copies of an element; they differ in what else is on that contig:
#   same_sequence_with_host          it also carries >= MIN_HOST_FLANK of unaligned sequence, so a
#                                    misjoined or rebinned contig that does include host sequence.
#   same_element_host_unassignable   >= MIN_ELEMENT_LEN with no host flank: a free-standing element
#                                    assembled in both samples. The element is real, but nothing
#                                    ties it to either genome, so the host call is unreliable
#                                    (this is the BMT104 vanB case).
#   undecidable_short_fragment       < MIN_ELEMENT_LEN with no host flank: too short to separate a
#                                    transferred element from the same contig binned twice
#                                    (the H37 case). Reported, never counted either way.
#   shared_element                   the ARG sits inside the longest block and BOTH contigs carry
#                                    >= MIN_HOST_FLANK outside it - two distinct genomes sharing
#                                    an element. This is the only positive transfer evidence here.
#   arg_outside_block                there is a shared block but neither ARG ORF is in it.
#   no_shared_block                  no alignment passing MIN_ALN / MIN_ID.
#
# Input : <4.17 output>/Donor_contig_queries.tsv (written by 4.17): one row per event x donor copy.
#         Uses Source_donor / Source_after to find each assembly, so no cohort is inferred from IDs.
# Output: Donor_contig_alignment.tsv next to the input - the query columns plus the measurements,
#         `verdict`, and `decidable`. 4.17 can gate its event list on it.
#
# Controls, which the log prints first: BMT104 vanB (51,173 bp donor vs 50,968 bp recipient) and
# H37 APH(3')-IIIa / SAT-4 (1.9 and 2.2 kb donor vs 10.8 kb recipient). Both are known rebinning
# artefacts and must NOT come out shared_element.
#
# Usage (submits itself; run from the login node):
#   bash donor_contig_alignment.sh <log_file> <path/to/Donor_contig_queries.tsv>

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: donor_contig_alignment.sh <log_file> <Donor_contig_queries.tsv>}
    q=$(realpath "${2:?Usage: donor_contig_alignment.sh <log_file> <Donor_contig_queries.tsv>}")
    sbatch -o "$log" --export=ALL,DC_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",DC_QUERIES="${q}" \
        "${SCRIPT_DIR}/donor_contig_alignment.sh"
    exit $?
fi
CONFIG_DIR=${DC_CONFIG_DIR:?}
QUERIES=${DC_QUERIES:?}
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
module load minimap2/2.28
echo "minimap2 $(minimap2 --version)"

# Source (as 4.17 writes it) -> datasetDir|evaluationDir
srcs=()
for c in Duke_short UNC_short Heston_short; do
    v=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$datasetDir|$evaluationDir")
    srcs+=("${c%_short}=$v")
done

work=${TMPDIR:-/tmp}/dc_${SLURM_JOB_ID}; mkdir -p "$work"; trap 'rm -rf "$work"' EXIT

$PY - "$QUERIES" "$work" "$THREADS" "${srcs[@]}" <<'PYEOF'
import csv, os, re, subprocess, sys
from collections import defaultdict
queries, work, threads = sys.argv[1:4]
src = {}
for a in sys.argv[4:]:
    k, v = a.split("=", 1)
    ds, ev = v.split("|")
    src[k] = {"ds": ds, "eval": ev}
out = os.path.join(os.path.dirname(queries), "Donor_contig_alignment.tsv")

MIN_ALN, MIN_ID = 300, 0.95        # alignment blocks kept
MIN_CONTAINED = 0.90               # shorter contig this aligned => one sequence, not two
MIN_HOST_FLANK = 2000              # unaligned sequence needed to call a flank real
MIN_ELEMENT_LEN = 20000            # a fully-matching contig this long is an element, not a fragment

def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

rows = read_tsv(queries)
print(f"{len(rows)} donor copies over {len({(r['Cohort'], r['PatientID'], r['ARO'], r['recipient_species']) for r in rows})} events")

# ---- contig sequences, read once per sample ----
need = defaultdict(set)
for r in rows:
    need[(r["Source_donor"], r["SampleID_before"])].add(r["contig_donor"])
    need[(r["Source_after"], r["SampleID_after"])].add(r["contig_after"])
seqs = {}
for (s, sid), want in sorted(need.items()):
    path = os.path.join(src[s]["eval"], f"{sid}_final_assembly.fasta")
    cur = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n = line[1:].split()[0]
                cur = n if n in want else None
                if cur:
                    seqs[(sid, cur)] = []
            elif cur:
                seqs[(sid, cur)].append(line.strip())
for k in list(seqs):
    seqs[k] = "".join(seqs[k])
missing = [(r["SampleID_before"], r["contig_donor"]) for r in rows if (r["SampleID_before"], r["contig_donor"]) not in seqs] \
        + [(r["SampleID_after"], r["contig_after"]) for r in rows if (r["SampleID_after"], r["contig_after"]) not in seqs]
if missing:
    sys.exit(f"{len(missing)} contigs not found, e.g. {missing[:3]}")

# ---- ORF coordinates, read once per sample ----
orfs = {}
for (s, sid), _ in sorted(need.items()):
    faa = os.path.join(src[s]["ds"], "5.5_AA_amr_assembly", sid, f"{sid}_genes.faa")
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                f = line[1:].split(" # ")
                orfs[(sid, f[0])] = (int(f[1]), int(f[2]))

def orf_span(sid, oid):
    v = orfs.get((sid, oid))
    return v if v else None

def merged(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out

def cov(iv, n):
    return sum(e - s for s, e in merged(iv)) / n if n else 0.0

res = []
for i, r in enumerate(rows, 1):
    dsid, dctg = r["SampleID_before"], r["contig_donor"]
    rsid, rctg = r["SampleID_after"], r["contig_after"]
    dseq, rseq = seqs[(dsid, dctg)], seqs[(rsid, rctg)]
    qf, tf = os.path.join(work, "q.fa"), os.path.join(work, "t.fa")
    with open(qf, "w") as fo:
        fo.write(f">donor\n{dseq}\n")
    with open(tf, "w") as fo:
        fo.write(f">recipient\n{rseq}\n")
    p = subprocess.run(["minimap2", "-x", "asm5", "-c", "-t", threads, tf, qf],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"minimap2 failed on row {i}: {p.stderr[-300:]}")
    hits = []
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if int(f[10]) >= MIN_ALN and int(f[9]) / int(f[10]) >= MIN_ID:
            hits.append({"qs": int(f[2]), "qe": int(f[3]), "ts": int(f[7]), "te": int(f[8]),
                         "id": int(f[9]) / int(f[10]), "blen": int(f[10])})
    dl, rl = len(dseq), len(rseq)
    short_len = min(dl, rl)
    o = {k: r[k] for k in rows[0]}
    d_orf, r_orf = orf_span(dsid, r["orf_donor"]), orf_span(rsid, r["orf_after"])
    if not hits:
        o.update({"donor_aln_frac": "0.00", "recipient_aln_frac": "0.00", "shorter_aln_frac": "0.00",
                  "best_identity": "", "longest_block_bp": 0, "arg_donor_in_block": "",
                  "arg_recipient_in_block": "", "donor_flank_bp": dl, "recipient_flank_bp": rl,
                  "shorter_flank_bp": short_len, "decidable": "TRUE", "verdict": "no_shared_block"})
    else:
        df = cov([(h["qs"], h["qe"]) for h in hits], dl)
        rf = cov([(h["ts"], h["te"]) for h in hits], rl)
        sf = df if dl <= rl else rf
        best = max(hits, key=lambda h: h["blen"])
        din = (d_orf and best["qs"] <= d_orf[0] - 1 and d_orf[1] <= best["qe"])
        rin = (r_orf and best["ts"] <= r_orf[0] - 1 and r_orf[1] <= best["te"])
        d_flank = dl - int(df * dl)
        r_flank = rl - int(rf * rl)
        s_flank = d_flank if dl <= rl else r_flank
        decidable = s_flank >= MIN_HOST_FLANK
        if sf >= MIN_CONTAINED:
            # The shorter contig is wholly contained in the longer, so the two samples hold one
            # sequence, not two genomes' copies. What that means depends on what else is on it:
            #   with host flank      -> a misjoined/rebinned contig that does carry host sequence
            #   long, no host flank  -> a free-standing element assembled in both samples; the
            #                           element is real but its host cannot be assigned (BMT104)
            #   short, no host flank -> a fragment; transfer and rebinning look identical
            if decidable:
                v = "same_sequence_with_host"
            elif short_len >= MIN_ELEMENT_LEN:
                v = "same_element_host_unassignable"
            else:
                v = "undecidable_short_fragment"
        elif (din or rin) and d_flank >= MIN_HOST_FLANK and r_flank >= MIN_HOST_FLANK:
            v = "shared_element"
        elif not (din or rin):
            v = "arg_outside_block"
        else:
            v = "ambiguous"
        o.update({"donor_aln_frac": f"{df:.2f}", "recipient_aln_frac": f"{rf:.2f}",
                  "shorter_aln_frac": f"{sf:.2f}", "best_identity": f"{best['id']:.4f}",
                  "longest_block_bp": best["blen"],
                  "arg_donor_in_block": "TRUE" if din else "FALSE",
                  "arg_recipient_in_block": "TRUE" if rin else "FALSE",
                  "donor_flank_bp": d_flank, "recipient_flank_bp": r_flank,
                  "shorter_flank_bp": s_flank,
                  "decidable": "TRUE" if decidable else "FALSE", "verdict": v})
    res.append(o)
    if i % 25 == 0 or i == len(rows):
        print(f"  {i}/{len(rows)} pairs", flush=True)

cols = list(rows[0]) + ["donor_aln_frac", "recipient_aln_frac", "shorter_aln_frac", "best_identity",
                        "longest_block_bp", "arg_donor_in_block", "arg_recipient_in_block",
                        "donor_flank_bp", "recipient_flank_bp", "shorter_flank_bp",
                        "decidable", "verdict"]
with open(out, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(res)

print("\n=== CONTROLS (both are known rebinning artefacts; neither may be shared_element) ===")
for o in res:
    if (o["PatientID"] == "BMT104" and "van" in o["Gene"]) or (o["PatientID"] == "H37"):
        print(f"  {o['PatientID']:<8}{o['Gene'][:22]:<23} donor {o['contig_length_donor']:>8}bp  "
              f"recipient {o['contig_length_after']:>8}bp  shorter_aln {o['shorter_aln_frac']:>5}  "
              f"flank {o['shorter_flank_bp']:>7}bp  -> {o['verdict']}")

print("\n=== verdicts, all rows ===")
cnt = defaultdict(int)
for o in res:
    cnt[o["verdict"]] += 1
for k in sorted(cnt, key=lambda x: -cnt[x]):
    print(f"  {k:<22} {cnt[k]}")
print("\n=== verdicts for the same-strain (both-tools) events only ===")
cnt2 = defaultdict(int)
for o in res:
    if o.get("strain") == "same strain":
        cnt2[o["verdict"]] += 1
for k in sorted(cnt2, key=lambda x: -cnt2[x]):
    print(f"  {k:<22} {cnt2[k]}")
se = [o for o in res if o["verdict"] == "shared_element"]
print(f"\n=== shared_element candidates ({len(se)}) ===")
for o in sorted(se, key=lambda x: -int(x["longest_block_bp"]))[:30]:
    print(f"  {o['Cohort']:<9} {o['PatientID']:<8} {o['recipient_species'][:24]:<25} {o['Gene'][:16]:<17} "
          f"from {o['family_donor'][:20]:<21} block {o['longest_block_bp']:>7}bp id {o['best_identity']} "
          f"flanks {o['donor_flank_bp']}/{o['recipient_flank_bp']} strain={o.get('strain','')}")
print(f"\nwrote {out}")
PYEOF
echo "done"
