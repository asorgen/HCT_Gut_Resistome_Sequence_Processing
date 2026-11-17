#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --job-name=tier2_mag_ani

# ANI between the "before" and "after" MAG of every tier-2 ARG gain event
# (bug-to-bug ARG transfer analysis, module 4.15 in HCT_Gut_Resistome_Species_Attribution).
#
# A tier-2 event is a species whose earlier MAG lacks a gene and whose later MAG has it.
# ANI separates a resident strain gaining the gene (ANI >= ~99.9%) from a new strain
# of the same species arriving with it (lower ANI, still >= 95% = same species).
#
# Input : Tier2_MAG_pairs_for_ANI.tsv (written by 4.15). One row per event (patient x
#         species x gene), so several rows can share a MAG pair. skani runs once per
#         unique pair, and the result is joined back to every row.
# MAGs  : <datasetDir>/3.3_reassemble_bins/<SampleID>/reassembled_bins/<MAG_version>.fa
#         datasetDir from each cohort's config. Cohort "Andermann" splits into
#         UNC_short (BMT* IDs) and Duke_short (the rest).
# Tool  : skani dist, ref = before MAG, query = after MAG. skani reports nothing below
#         ~80% ANI; those pairs get ANI = NA and skani_note = below_skani_cutoff.
# Output: Tier2_MAG_pairs_ANI.tsv next to the input: every input column plus
#         ANI, AF_before (aligned fraction of the before MAG), AF_after, skani_note.
#
# Usage (submits itself; run from the login node):
#   bash tier2_mag_ani.sh <log_file> <path/to/Tier2_MAG_pairs_for_ANI.tsv>

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: tier2_mag_ani.sh <log_file> <pairs_tsv>}
    pairs=$(realpath "${2:?Usage: tier2_mag_ani.sh <log_file> <pairs_tsv>}")
    sbatch -o "$log" --export=ALL,ANI_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",ANI_PAIRS="${pairs}" \
        "${SCRIPT_DIR}/tier2_mag_ani.sh"
    exit $?
fi
CONFIG_DIR=${ANI_CONFIG_DIR:?}
PAIRS=${ANI_PAIRS:?}
source "${CONFIG_DIR}/private.config"
# Full path, not `conda activate`: env names collide between miniforge3 and .conda on this HPC
SKANI=${HPC_HOME}/.conda/envs/bin_refine/bin/skani
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
OUT=$(dirname "$PAIRS")/Tier2_MAG_pairs_ANI.tsv

dirs=()
for c in Duke_short UNC_short Heston_short; do
    d=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$datasetDir")
    [[ -d $d ]] || { echo "No datasetDir for $c ($d)"; exit 1; }
    dirs+=("${c%_short}=$d")
done
echo "skani: $($SKANI --version)"
echo "pairs: $PAIRS"

$PY - "$SKANI" "$THREADS" "$PAIRS" "$OUT" "${dirs[@]}" <<'PYEOF'
import csv, os, subprocess, sys, tempfile
skani, threads, pairs_tsv, out = sys.argv[1:5]
dirs = dict(a.split("=", 1) for a in sys.argv[5:])

def cohort(row, sid):
    if row["Cohort"] == "Heston":
        return "Heston"
    return "UNC" if sid.startswith("BMT") else "Duke"

def mag_path(row, w):
    sid = row[f"SampleID_{w}"]
    return os.path.join(dirs[cohort(row, sid)], "3.3_reassemble_bins", sid,
                        "reassembled_bins", row[f"MAG_version_{w}"] + ".fa")

with open(pairs_tsv, newline="") as fh:
    rdr = csv.DictReader(fh, delimiter="\t")
    cols = rdr.fieldnames
    rows = list(rdr)

missing = sorted({p for r in rows for p in (mag_path(r, "before"), mag_path(r, "after"))
                  if not os.path.exists(p)})
if missing:
    sys.exit(f"{len(missing)} MAG fasta(s) missing, e.g. {missing[:3]}")

# accessory events first (the TaskNote priority); skani each unique pair once
rows_sorted = sorted(rows, key=lambda r: r.get("tier2_accessory") != "TRUE")
res = {}
with tempfile.TemporaryDirectory() as tmp:
    for r in rows_sorted:
        key = (mag_path(r, "before"), mag_path(r, "after"))
        if key in res:
            continue
        o = os.path.join(tmp, "o.tsv")
        p = subprocess.run([skani, "dist", "-t", threads, "-q", key[1], "-r", key[0], "-o", o],
                           capture_output=True, text=True)
        if p.returncode != 0:
            res[key] = ("NA", "NA", "NA", "skani_error:" + p.stderr.strip().splitlines()[-1][:80])
            continue
        hits = list(csv.DictReader(open(o), delimiter="\t"))
        if not hits:
            res[key] = ("NA", "NA", "NA", "below_skani_cutoff")
        else:
            h = hits[0]
            res[key] = (h["ANI"], h["Align_fraction_ref"], h["Align_fraction_query"], "")
        print(f"  {len(res):>4}  {r['Cohort']:<9} {r['gtdb_species']:<35} ANI {res[key][0]:>6}  "
              f"AF {res[key][1]}/{res[key][2]}  {res[key][3]}", flush=True)

with open(out, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(cols + ["ANI", "AF_before", "AF_after", "skani_note"])
    for r in rows:   # input order, so the file lines up with 4.15's
        w.writerow([r[c] for c in cols] + list(res[(mag_path(r, "before"), mag_path(r, "after"))]))

ani = [float(v[0]) for v in res.values() if v[0] != "NA"]
print(f"\n{out}: {len(rows)} rows, {len(res)} unique MAG pairs, {len(ani)} with ANI")
for lab, lo, hi in ((">= 99.9", 99.9, 101), ("99.0-99.9", 99.0, 99.9), ("95-99", 95.0, 99.0), ("< 95", 0, 95.0)):
    print(f"  ANI {lab:<10} {sum(lo <= a < hi for a in ani)}")
PYEOF
echo "done"
