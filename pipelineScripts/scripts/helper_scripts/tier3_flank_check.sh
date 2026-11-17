#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=02:00:00
#SBATCH --job-name=tier3_flank_check

# Insertion / flank check for tier-3 ARG gains (bug-to-bug ARG transfer, module 4.15).
#
# A tier-3 event is a same-strain gain (skani ANI >= 99.9% between the before and after
# MAG). The most direct evidence that the resident chromosome acquired the gene is an
# insertion: the sequence flanking the ARG region on the "after" gene contig is present
# in the "before" genome, while the ARG region itself is absent, and the two flanks sit
# next to each other there.
#
# Events checked (from Tier2_events.tsv):
#   * tier3 & host_anchored (gene contig >= 25 kb and Kraken agrees with the MAG)
#   * tier3 & gene in the aadE-sat4-aphA-3 cluster (aad(6), APH(3')-IIIa, SAT-4)
#   * element_shared_with set (one element binned to two species, e.g. BMT104 vanB),
#     tier3 or not, plus the same patient's other-species rows for that gene
#
# For each event, the after gene contig is aligned (minimap2 -x asm5) to
#   (a) the before MAG (3.3 reassembled bin, MAG_version from the MAG catalog), and
#   (b) the before sample's whole assembly (1.2), which catches flanks left unbinned.
# Gene span = the event gene's ORF(s) on the contig (RGI hit >= 80% identity, Prodigal
# coordinates from 5.5 <ID>_genes.faa). Alignments kept: >= 300 bp, >= 95% identity.
# Inserted segment = the unaligned stretch of the contig around the gene, bounded by
# the nearest aligned blocks. A flank is present if an aligned block bounds the segment
# on that side; if the segment runs to the contig end, FALSE + note contig_end_*.
# clean_junction = the two bounding blocks hit the same target contig, same strand,
# with their facing ends <= 1 kb apart (the flanks are adjacent in "before").
# gene_*_before_asm = a sensitive gene-only search (minimap2 asm20) of the whole before
# assembly with the RGI call on the hit contig: catches the gene under another RGI name
# (sister allele, e.g. CfxA2 vs CfxA3), which tier 2's exact-ARO absence check misses.
# Mobile-element genes in the inserted segment (or the whole contig if no flank aligns)
# come from the existing 5.5 bakta protein annotation of the after sample.
#
# Output (next to Tier2_events.tsv):
#   Tier3_insertion_check.tsv   one row per event: keys + before_flank_left_present,
#                               before_flank_right_present, gene_region_in_before,
#                               inserted_len, notes; the same for the whole before
#                               assembly (_asm columns); mobile genes in the segment
#   Tier3_insertion_check_BMT104.tsv   the two BMT104 vanB contigs against each other
#                               and against the BMT104 C. leptum / Ruthenibacterium MAGs
#
# Usage (submits itself; run from the login node):
#   bash tier3_flank_check.sh <log_file> <path/to/Tier2_events.tsv>

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: tier3_flank_check.sh <log_file> <Tier2_events.tsv>}
    events=$(realpath "${2:?Usage: tier3_flank_check.sh <log_file> <Tier2_events.tsv>}")
    sbatch -o "$log" --export=ALL,FC_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",FC_EVENTS="${events}" \
        "${SCRIPT_DIR}/tier3_flank_check.sh"
    exit $?
fi
CONFIG_DIR=${FC_CONFIG_DIR:?}
EVENTS=${FC_EVENTS:?}
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
module load minimap2/2.28
echo "minimap2 $(minimap2 --version)"

# cohort=datasetDir|evaluationDir|PROCESSED_ROOT/<dataset>_tables
cohorts=()
for c in Duke_short UNC_short Heston_short; do
    v=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$datasetDir|$evaluationDir|$PROCESSED_ROOT/${dataset}_tables/${dataset}")
    cohorts+=("${c%_short}=$v")
done

work=${TMPDIR:-/tmp}/fc_${SLURM_JOB_ID}
mkdir -p "$work"; trap 'rm -rf "$work"' EXIT

$PY - "$EVENTS" "$work" "$THREADS" "${cohorts[@]}" <<'PYEOF'
import csv, os, re, subprocess, sys
from collections import defaultdict
events_tsv, work, threads = sys.argv[1:4]
coh = {}
for a in sys.argv[4:]:
    k, v = a.split("=", 1)
    ds, ev, tab = v.split("|")
    coh[k] = {"ds": ds, "eval": ev, "tab": tab}
outdir = os.path.dirname(events_tsv)
CLUSTER = {"aad(6)", "APH(3')-IIIa", "SAT-4"}
MOBILE = re.compile(r"transposase|integrase|recombinase|resolvase|invertase|excisionase|relaxase|"
                    r"\bMob[A-Z]?\b|mobili[sz]ation|conjugal|conjugative|Tn916|Tn5397|Tn1549|"
                    r"\bIS\d|insertion sequence|phage|plasmid", re.I)
FLANK, MIN_WIN, MIN_ALN, MIN_ID, ADJ = 5000, 500, 300, 0.95, 1000

def cohort_of(row, sid):
    return "Heston" if row["Cohort"] == "Heston" else ("UNC" if sid.startswith("BMT") else "Duke")

def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

# ---- events to check ----
ev = read_tsv(events_tsv)
shared = {(r["PatientID"], r["Gene"]) for r in ev if r.get("element_shared_with")}
sel = [r for r in ev if (r["tier3"] == "TRUE" and (r["host_anchored"] == "TRUE" or r["Gene"] in CLUSTER))
       or ((r["PatientID"], r["Gene"]) in shared and r["tier2_accessory"] == "TRUE")]
print(f"{len(sel)} event rows selected")

# ---- MAG versions ----
magv = {}
for c in coh:
    for r in read_tsv(coh[c]["tab"] + "_mag_catalog.tsv"):
        magv[(c, r["SampleID"], r["MAG"])] = r["MAG_version"]

def mag_fa(c, sid, mag):
    return os.path.join(coh[c]["ds"], "3.3_reassemble_bins", sid, "reassembled_bins", magv[(c, sid, mag)] + ".fa")

def asm_fa(c, sid):
    return os.path.join(coh[c]["eval"], f"{sid}_final_assembly.fasta")

# ---- contig sequence, ORF coordinates, RGI hits, bakta products for an after sample ----
_cache = {}
def sample_info(c, sid, contig):
    key = (c, sid, contig)
    if key in _cache:
        return _cache[key]
    seq, keep = [], False
    with open(asm_fa(c, sid)) as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = line[1:].split()[0] == contig
            elif keep:
                seq.append(line.strip())
    seq = "".join(seq)
    if not seq:
        sys.exit(f"contig {contig} not found in {asm_fa(c, sid)}")
    d55 = os.path.join(coh[c]["ds"], "5.5_AA_amr_assembly", sid)
    orfs = {}   # orf_id -> (start, end)
    with open(os.path.join(d55, f"{sid}_genes.faa")) as fh:
        for line in fh:
            if line.startswith(">" + contig + "_"):
                f = line[1:].split(" # ")
                orfs[f[0]] = (int(f[1]), int(f[2]))
    rgi = defaultdict(list)  # Best_Hit_ARO -> [(start, end, identity)]
    for r in read_tsv(os.path.join(d55, f"{sid}.rgi.txt")):
        oid = r["ORF_ID"].split(" # ")[0]
        if oid in orfs:
            rgi[r["Best_Hit_ARO"]].append((*orfs[oid], float(r["Best_Identities"])))
    prod = {}
    bt = os.path.join(d55, "bakta", f"{sid}.tsv")
    if os.path.exists(bt):
        with open(bt) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[0] in orfs and len(f) > 3:
                    prod[f[0]] = f[3]
    _cache[key] = (seq, orfs, rgi, prod)
    return _cache[key]

# ---- minimap2 ----
_aln = {}
def align(qname, qseq, target):
    key = (qname, target)
    if key in _aln:
        return _aln[key]
    q = os.path.join(work, "q.fa")
    with open(q, "w") as fo:
        fo.write(f">{qname}\n{qseq}\n")
    p = subprocess.run(["minimap2", "-x", "asm5", "-c", "-t", threads, target, q],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"minimap2 failed on {target}: {p.stderr[-300:]}")
    hits = []
    for line in p.stdout.splitlines():
        f = line.split("\t")
        qs, qe, tname, tlen, ts, te, match, blen = int(f[2]), int(f[3]), f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if blen >= MIN_ALN and match / blen >= MIN_ID:
            hits.append({"qs": qs, "qe": qe, "t": tname, "tlen": tlen, "ts": ts, "te": te,
                         "strand": f[4], "id": match / blen})
    _aln[key] = hits
    return hits

def covered(hits, a, b):
    """bp of [a, b) covered by the union of query intervals."""
    if b <= a:
        return 0
    iv = sorted((max(h["qs"], a), min(h["qe"], b)) for h in hits if h["qe"] > a and h["qs"] < b)
    tot, cur_s, cur_e = 0, None, None
    for s, e in iv:
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot

def flank_metrics(hits, L, gs, ge):
    """gs/ge 0-based half-open gene span on a contig of length L.

    The inserted segment is the unaligned stretch around the gene, bounded by the
    nearest aligned blocks. A flank is the sequence beyond the segment's edge (not
    beyond the gene's: the element is usually longer than the gene). Flank present =
    an aligned block bounds the segment on that side; if the segment runs to the
    contig end there is no flank sequence on the contig to test (FALSE, noted
    contig_end_left/right). flank_*_frac = aligned share of up to 5 kb beyond the edge.
    """
    gfrac = covered(hits, gs, ge) / (ge - gs)
    left = [h for h in hits if h["qe"] <= gs + 50]
    right = [h for h in hits if h["qs"] >= ge - 50]
    lb = max(left, key=lambda h: h["qe"]) if left else None
    rb = min(right, key=lambda h: h["qs"]) if right else None
    ins_s = lb["qe"] if lb else 0
    ins_e = rb["qs"] if rb else L
    def frac(a, b):
        return f"{covered(hits, a, b) / (b - a):.2f}" if b - a >= MIN_WIN else ""
    lp, lf = ("TRUE" if lb else "FALSE"), frac(max(0, ins_s - FLANK), ins_s)
    rp, rf = ("TRUE" if rb else "FALSE"), frac(ins_e, min(L, ins_e + FLANK))
    notes = []
    if not lb: notes.append("contig_end_left")
    if not rb: notes.append("contig_end_right")
    junction, tgap, tname = "", "", ""
    if lb and rb:
        tname = lb["t"] if lb["t"] == rb["t"] else f"{lb['t']}/{rb['t']}"
        if lb["t"] == rb["t"]:
            # the flank ends that face the insertion, in target coordinates
            lend = lb["te"] if lb["strand"] == "+" else lb["ts"]
            rend = rb["ts"] if rb["strand"] == "+" else rb["te"]
            tgap = abs(rend - lend)
            junction = "TRUE" if (tgap <= ADJ and lb["strand"] == rb["strand"]) else "FALSE"
        else:
            junction = "FALSE"
            notes.append("flanks_on_different_target_contigs")
    return {"flank_left_present": lp, "flank_left_frac": lf, "flank_right_present": rp,
            "flank_right_frac": rf, "gene_region_frac": f"{gfrac:.2f}",
            "gene_region_in_before": "TRUE" if gfrac >= 0.5 else "FALSE",
            "inserted_start": ins_s, "inserted_end": ins_e, "inserted_len": ins_e - ins_s,
            "clean_junction": junction, "junction_target_gap": tgap, "junction_target": tname,
            "contig_aligned_frac": f"{covered(hits, 0, L) / L:.2f}", "notes": ";".join(notes)}

_rgi_before = {}
def gene_in_before(c, sid, gseq):
    """Sensitive gene-only search of the whole before assembly (minimap2 asm20), to catch
    the gene under a different RGI name (sister allele, mosaic) that tier 2's exact-ARO
    absence check misses. Returns coverage of the gene, identity, hit contig, RGI there."""
    q = os.path.join(work, "g.fa")
    with open(q, "w") as fo:
        fo.write(f">g\n{gseq}\n")
    p = subprocess.run(["minimap2", "-x", "asm20", "-c", asm_fa(c, sid), q], capture_output=True, text=True)
    hits = [l.split("\t") for l in p.stdout.splitlines()]
    if not hits:
        return {"gene_cov_before_asm": "0.00", "gene_id_before_asm": "", "gene_hit_contig_before": "",
                "gene_hit_rgi_before": ""}
    if (c, sid) not in _rgi_before:
        d = defaultdict(list)
        for r in read_tsv(os.path.join(coh[c]["ds"], "5.5_AA_amr_assembly", sid, f"{sid}.rgi.txt")):
            d[r["ORF_ID"].split(" # ")[0].rsplit("_", 1)[0]].append(f"{r['Best_Hit_ARO']}({r['Best_Identities']})")
        _rgi_before[(c, sid)] = d
    hd = [{"qs": int(h[2]), "qe": int(h[3])} for h in hits]
    best = max(hits, key=lambda h: int(h[10]))
    tg = sorted({h[5] for h in hits})
    return {"gene_cov_before_asm": f"{covered(hd, 0, len(gseq)) / len(gseq):.2f}",
            "gene_id_before_asm": f"{int(best[9]) / int(best[10]):.3f}",
            "gene_hit_contig_before": ",".join(tg)[:300],
            "gene_hit_rgi_before": "; ".join(", ".join(_rgi_before[(c, sid)].get(t, ["none"])) for t in tg)[:300]}

def mobile_in(orfs, prod, a, b):
    hits = sorted({prod[o] for o, (s, e) in orfs.items() if o in prod and s - 1 >= a and e <= b
                   and MOBILE.search(prod[o])})
    return len(hits), " | ".join(hits)[:400]

# ---- per event ----
keys = ["Cohort", "PatientID", "Gene", "ARO", "gtdb_species", "tier3", "host_anchored", "ANI",
        "SampleID_before", "MAG_before", "SampleID_after", "MAG_after", "contig_after", "contig_length_after"]
mcols = ["flank_left_present", "flank_left_frac", "flank_right_present", "flank_right_frac",
         "gene_region_frac", "gene_region_in_before", "inserted_start", "inserted_end", "inserted_len",
         "clean_junction", "junction_target_gap", "junction_target", "contig_aligned_frac", "notes"]
out_cols = keys + ["gene_start", "gene_end", "gene_identity"] + \
           ["before_" + c for c in mcols] + [c + "_asm" for c in mcols] + \
           ["gene_cov_before_asm", "gene_id_before_asm", "gene_hit_contig_before", "gene_hit_rgi_before",
            "n_mobile_genes_in_segment", "mobile_genes_in_segment"]
rows_out = []
for r in sel:
    ca, cb = cohort_of(r, r["SampleID_after"]), cohort_of(r, r["SampleID_before"])
    contig = r["contig_after"]
    seq, orfs, rgi, prod = sample_info(ca, r["SampleID_after"], contig)
    g = [h for h in rgi.get(r["Gene"], []) if h[2] >= 80]
    if not g:
        print(f"  WARN no >=80% RGI hit for {r['Gene']} on {r['SampleID_after']} {contig}")
        continue
    gs, ge = min(h[0] for h in g) - 1, max(h[1] for h in g)
    qname = f"{r['SampleID_after']}|{contig}"
    m_mag = flank_metrics(align(qname, seq, mag_fa(cb, r["SampleID_before"], r["MAG_before"])), len(seq), gs, ge)
    m_asm = flank_metrics(align(qname, seq, asm_fa(cb, r["SampleID_before"])), len(seq), gs, ge)
    # segment for the mobile-gene scan: from the whole before assembly (the fuller picture)
    nmob, mob = mobile_in(orfs, prod, m_asm["inserted_start"], m_asm["inserted_end"])
    gb = gene_in_before(cb, r["SampleID_before"], seq[gs:ge])
    o = {k: r[k] for k in keys}
    o.update({"gene_start": gs + 1, "gene_end": ge, "gene_identity": max(h[2] for h in g)})
    o.update({"before_" + k: v for k, v in m_mag.items()})
    o.update({k + "_asm": v for k, v in m_asm.items()})
    o.update(gb)
    o.update({"n_mobile_genes_in_segment": nmob, "mobile_genes_in_segment": mob})
    rows_out.append(o)
    print(f"  {r['PatientID']:<7} {r['Gene'][:22]:<22} {r['gtdb_species'][:26]:<26} "
          f"MAG: L={m_mag['flank_left_present']:<5} R={m_mag['flank_right_present']:<5} "
          f"| asm: L={m_asm['flank_left_present']:<5} R={m_asm['flank_right_present']:<5} "
          f"ins={m_asm['inserted_len']} junction={m_asm['clean_junction'] or 'NA'} "
          f"| gene in before asm: cov={gb['gene_cov_before_asm']} {gb['gene_hit_rgi_before'][:40]} | mobile={nmob}", flush=True)

with open(os.path.join(outdir, "Tier3_insertion_check.tsv"), "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=out_cols, delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(rows_out)
print(f"wrote Tier3_insertion_check.tsv: {len(rows_out)} rows")

# ---- BMT104: the two vanB contigs vs each other and vs the C. leptum / Ruthenibacterium MAGs ----
bmt = [("BMT104D24", "NODE_374_length_51173_cov_27.573262"), ("BMT104D29", "NODE_338_length_50968_cov_25.521930")]
targets = [("BMT104D-53", "bin.15", "C. leptum, before"), ("BMT104D24", "bin.2", "C. leptum, day 24 (holds D24 contig)"),
           ("BMT104D14", "bin.4", "Ruthenibacterium, day 14 (before)"),
           ("BMT104D24", "bin.21", "Ruthenibacterium, day 24 (lacks vanB)"),
           ("BMT104D29", "bin.14", "Ruthenibacterium, day 29 (holds D29 contig)")]
brow = []
for sid, contig in bmt:
    seq, orfs, rgi, prod = sample_info("UNC", sid, contig)
    g = [h for k, v in rgi.items() if k == "vanB" or "vanB cluster" in k for h in v if h[2] >= 80]
    gs, ge = min(h[0] for h in g) - 1, max(h[1] for h in g)
    qname = f"{sid}|{contig}"
    tl = [(t_sid, mag, lab, mag_fa("UNC", t_sid, mag)) for t_sid, mag, lab in targets]
    other = [x for x in bmt if x[0] != sid][0]
    ofa = os.path.join(work, f"{other[0]}_vanB_contig.fa")
    with open(ofa, "w") as fo:
        fo.write(f">{other[0]}|{other[1]}\n{sample_info('UNC', *other)[0]}\n")
    tl.append((other[0], other[1], "the other BMT104 vanB contig", ofa))
    for t_sid, mag, lab, fa in tl:
        hits = align(qname, seq, fa)
        m = flank_metrics(hits, len(seq), gs, ge)
        brow.append({"query_sample": sid, "query_contig": contig, "van_operon_start": gs + 1, "van_operon_end": ge,
                     "target_sample": t_sid, "target": mag, "target_label": lab,
                     "mean_identity": f"{sum(h['id'] for h in hits) / len(hits):.4f}" if hits else "", **m})
        print(f"  BMT104 {sid} contig vs {t_sid} {mag:<36} aligned={m['contig_aligned_frac']} "
              f"L={m['flank_left_present']} R={m['flank_right_present']} van={m['gene_region_frac']} "
              f"junction={m['clean_junction'] or 'NA'} ({lab})", flush=True)
with open(os.path.join(outdir, "Tier3_insertion_check_BMT104.tsv"), "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=list(brow[0]), delimiter="\t", lineterminator="\n")
    w.writeheader(); w.writerows(brow)
print("wrote Tier3_insertion_check_BMT104.tsv")
PYEOF
echo "done"
