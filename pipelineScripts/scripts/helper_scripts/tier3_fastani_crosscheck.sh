#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --job-name=tier3_fastani

# fastANI cross-check of the skani ANI behind tier 3 (bug-to-bug ARG transfer, module 4.15).
#
# Tier 3 turns on ANI >= 99.9% between the before and after MAG, but 99.9% is close to the limit
# of what either tool resolves on draft MAGs, so the call should not rest on one estimator. This
# runs fastANI over the same MAG pairs and reports both numbers side by side, flagging any pair
# where the two tools fall on opposite sides of the threshold.
#
# fastANI is a different method (fixed 3 kb fragments, MashMap mapping) from skani (sparse k-mer
# sketching), so agreement is meaningful and disagreement localises the uncertainty.
#
# By default it checks only the pairs that matter for the call: tier-3 pairs and anything in the
# unresolved band. --all does every pair in the file.
#
# Input : Tier2_MAG_pairs_ANI.tsv (skani results, from tier2_mag_ani.sh)
#         Tier2_events.tsv for the tier3 flag (optional; without it every pair is checked)
# Output: Tier3_fastANI_crosscheck.tsv next to the input - the pair keys, skani ANI and aligned
#         fractions, fastANI ANI, its orthologous-fragment counts, the difference, and agree/differ
#         against the 99.9% threshold.
#
# Usage (submits itself; run from the login node):
#   bash tier3_fastani_crosscheck.sh <log_file> <path/to/Tier2_MAG_pairs_ANI.tsv> [--all]

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: tier3_fastani_crosscheck.sh <log_file> <Tier2_MAG_pairs_ANI.tsv> [--all]}
    ani=$(realpath "${2:?Usage: tier3_fastani_crosscheck.sh <log_file> <Tier2_MAG_pairs_ANI.tsv> [--all]}")
    sbatch -o "$log" --export=ALL,FA_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",FA_ANI="${ani}",FA_ALL="${3:-}" \
        "${SCRIPT_DIR}/tier3_fastani_crosscheck.sh"
    exit $?
fi
CONFIG_DIR=${FA_CONFIG_DIR:?}
ANI_FILE=${FA_ANI:?}
ALL=${FA_ALL:-}
source "${CONFIG_DIR}/private.config"
# Full path, not `conda activate`: env names collide between miniforge3 and .conda on this HPC
FASTANI=${HPC_HOME}/.conda/envs/bin_refine/bin/fastANI
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
echo "fastANI $($FASTANI --version 2>&1 | head -1)"

dirs=()
for c in Duke_short UNC_short Heston_short; do
    v=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; \
        echo "$datasetDir|$PROCESSED_ROOT/${dataset}_tables/${dataset}_mag_catalog.tsv")
    dirs+=("${c%_short}=$v")
done

work=${TMPDIR:-/tmp}/fa_${SLURM_JOB_ID}; mkdir -p "$work"; trap 'rm -rf "$work"' EXIT

$PY - "$FASTANI" "$THREADS" "$ANI_FILE" "$work" "$ALL" "${dirs[@]}" <<'PYEOF'
import csv, os, subprocess, sys
fastani, threads, ani_file, work, allflag = sys.argv[1:6]
coh = {}
for a in sys.argv[6:]:
    k, v = a.split("=", 1)
    ds, cat = v.split("|")
    coh[k] = {"ds": ds, "cat": cat}
dirs = {k: v["ds"] for k, v in coh.items()}
outdir = os.path.dirname(ani_file)
SAME = 99.9

def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def cohort_of(row, sid):
    return "Heston" if row["Cohort"] == "Heston" else ("UNC" if sid.startswith("BMT") else "Duke")

rows = read_tsv(ani_file)
# MAG_version_* is not in the ANI file, so take it from the MAG catalogs, as tier2_mag_ani.sh does
magv = {}
for c, v in coh.items():
    if not os.path.exists(v["cat"]):
        sys.exit(f"missing MAG catalog for {c}: {v['cat']}")
    for r in read_tsv(v["cat"]):
        magv[(c, r["SampleID"], r["MAG"])] = r["MAG_version"]

def mag_fa(c, sid, mag):
    v = magv.get((c, sid, mag))
    if v is None:
        return None
    return os.path.join(dirs[c], "3.3_reassemble_bins", sid, "reassembled_bins", v + ".fa")

# which pairs to check
events = {}
ev_path = os.path.join(outdir, "Tier2_events.tsv")
if os.path.exists(ev_path):
    for e in read_tsv(ev_path):
        k = (e["Cohort"], e["PatientID"], e["ARO"], e["gtdb_species"],
             e["SampleID_before"], e["MAG_before"], e["SampleID_after"], e["MAG_after"])
        events[k] = e.get("tier3", "")

def wanted(r):
    if allflag == "--all":
        return True
    try:
        a = float(r["ANI"])
    except (ValueError, KeyError):
        return True
    k = (r["Cohort"], r["PatientID"], r["ARO"], r["gtdb_species"],
         r["SampleID_before"], r["MAG_before"], r["SampleID_after"], r["MAG_after"])
    return events.get(k) == "TRUE" or a >= 99.0     # tier 3, or in / above the unresolved band

pairs, meta = {}, {}
for r in rows:
    if not wanted(r):
        continue
    cb, ca = cohort_of(r, r["SampleID_before"]), cohort_of(r, r["SampleID_after"])
    fb, fa = mag_fa(cb, r["SampleID_before"], r["MAG_before"]), mag_fa(ca, r["SampleID_after"], r["MAG_after"])
    if not fb or not fa or not os.path.exists(fb) or not os.path.exists(fa):
        continue
    key = (fb, fa)
    pairs.setdefault(key, r)
    meta.setdefault(key, []).append(r)
print(f"{len(rows)} rows in the skani file; {len(pairs)} unique MAG pairs to cross-check"
      f"{' (all)' if allflag == '--all' else ' (tier 3 or ANI >= 99)'}")

# fastANI in one batch: --ql / --rl line up pairwise only if run per pair, so run per pair but
# reuse the same process pool size; each call is ~1 s on MAG-sized genomes.
res = {}
for n, (fb, fa) in enumerate(pairs, 1):
    o = os.path.join(work, "o.txt")
    p = subprocess.run([fastani, "-q", fa, "-r", fb, "-o", o, "-t", threads],
                       capture_output=True, text=True)
    if p.returncode != 0:
        res[(fb, fa)] = ("", "", "", "fastani_error")
        continue
    line = open(o).read().split("\n")[0].split("\t") if os.path.getsize(o) else []
    if len(line) >= 5:
        res[(fb, fa)] = (f"{float(line[2]):.4f}", line[3], line[4], "")
    else:
        res[(fb, fa)] = ("", "", "", "below_fastani_cutoff")
    if n % 20 == 0 or n == len(pairs):
        print(f"  {n}/{len(pairs)} pairs", flush=True)

keys = ["Cohort", "PatientID", "Gene", "ARO", "gtdb_species", "SampleID_before", "MAG_before",
        "SampleID_after", "MAG_after"]
out = os.path.join(outdir, "Tier3_fastANI_crosscheck.tsv")
cols = keys + ["skani_ANI", "AF_before", "AF_after", "fastANI", "fastANI_frags_mapped",
               "fastANI_frags_total", "delta_fastANI_minus_skani", "verdict", "note"]
n_out = 0
disagree = []
with open(out, "w", newline="") as fo:
    w = csv.DictWriter(fo, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for key, rs in meta.items():
        f_ani, fm, ft, note = res.get(key, ("", "", "", "not_run"))
        for r in rs:
            try:
                s = float(r["ANI"])
            except ValueError:
                s = None
            d = f"{float(f_ani) - s:+.4f}" if (f_ani and s is not None) else ""
            if f_ani and s is not None:
                v = "agree" if (float(f_ani) >= SAME) == (s >= SAME) else "DIFFER"
            else:
                v = ""
            if v == "DIFFER":
                disagree.append((r, s, float(f_ani)))
            w.writerow({**{k: r.get(k, "") for k in keys},
                        "skani_ANI": r.get("ANI", ""), "AF_before": r.get("AF_before", ""),
                        "AF_after": r.get("AF_after", ""), "fastANI": f_ani,
                        "fastANI_frags_mapped": fm, "fastANI_frags_total": ft,
                        "delta_fastANI_minus_skani": d, "verdict": v, "note": note})
            n_out += 1
print(f"\n{out}: {n_out} rows over {len(pairs)} MAG pairs")
ok = [(k, v) for k, v in res.items() if v[0]]
print(f"  fastANI returned a value for {len(ok)} of {len(pairs)} pairs")
if ok:
    ds = [float(v[0]) - float(pairs[k]["ANI"]) for k, v in ok if pairs[k]["ANI"]]
    print(f"  fastANI minus skani: median {sorted(ds)[len(ds)//2]:+.4f}, "
          f"range {min(ds):+.4f} to {max(ds):+.4f}")
print(f"  pairs where the two tools land on opposite sides of {SAME}%: {len({id(r) for r, _, _ in disagree})} rows")
for r, s, f in disagree[:20]:
    print(f"    {r['Cohort']:<9} {r['PatientID']:<8} {r['gtdb_species'][:26]:<26} {r['Gene'][:16]:<16} "
          f"skani {s:7.3f}  fastANI {f:7.3f}")
PYEOF
echo "done"
