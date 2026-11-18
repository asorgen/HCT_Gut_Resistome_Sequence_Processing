#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=02:00:00
#SBATCH --job-name=tier2_mag_fastani

# fastANI cross-check of the skani strain calls for tier-2 ARG gain events
# (bug-to-bug ARG transfer analysis, module 4.15 in HCT_Gut_Resistome_Species_Attribution).
#
# tier2_mag_ani.sh (skani) calls a gain "same strain" at ANI >= 99.9%, which is close to
# skani's precision on draft MAGs. This repeats the before / after MAG comparison with fastANI,
# an independent method, for every accessory event (tier2_accessory = TRUE), so the flagship
# calls (D20248, H04, H75, H61) and the 99-99.9% "unresolved" group can be checked against it.
#
# Input : Tier2_MAG_pairs_for_ANI.tsv (written by 4.15) and Tier2_MAG_pairs_ANI.tsv (skani,
#         from tier2_mag_ani.sh), both in the 4.15 output directory. Only rows with
#         tier2_accessory = TRUE are run; fastANI runs once per unique MAG pair.
# MAGs  : <datasetDir>/3.3_reassemble_bins/<SampleID>/reassembled_bins/<MAG_version>.fa,
#         as in tier2_mag_ani.sh.
# Tool  : fastANI 1.32 (bin_refine env), both directions (before -> after and after -> before),
#         because fastANI is asymmetric. fastANI reports nothing below ~80% ANI.
# Output: Tier2_MAG_pairs_fastANI.tsv next to the input: the event keys, the skani ANI,
#         fastANI ANI in each direction and their mean, matched / total fragments in each
#         direction, and fastani_note. Then a skani-vs-fastANI table of the strain calls.
#
# Usage (submits itself; run from the login node):
#   bash tier2_mag_fastani.sh <log_file> <path/to/Tier2_MAG_pairs_for_ANI.tsv>

set -uo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    log=${1:?Usage: tier2_mag_fastani.sh <log_file> <pairs_tsv>}
    pairs=$(realpath "${2:?Usage: tier2_mag_fastani.sh <log_file> <pairs_tsv>}")
    sbatch -o "$log" --export=ALL,ANI_CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)",ANI_PAIRS="${pairs}" \
        "${SCRIPT_DIR}/tier2_mag_fastani.sh"
    exit $?
fi
CONFIG_DIR=${ANI_CONFIG_DIR:?}
PAIRS=${ANI_PAIRS:?}
source "${CONFIG_DIR}/private.config"
# Full path, not `conda activate`: env names collide between miniforge3 and .conda on this HPC
FASTANI=${HPC_HOME}/.conda/envs/bin_refine/bin/fastANI
PY=/usr/bin/python3
THREADS=${SLURM_CPUS_PER_TASK:-4}
SKANI_OUT=$(dirname "$PAIRS")/Tier2_MAG_pairs_ANI.tsv
OUT=$(dirname "$PAIRS")/Tier2_MAG_pairs_fastANI.tsv
[[ -s $SKANI_OUT ]] || { echo "No skani results at $SKANI_OUT: run tier2_mag_ani.sh first"; exit 1; }

dirs=()
for c in Duke_short UNC_short Heston_short; do
    d=$(source ${CONFIG_DIR}/${c}-read.config >/dev/null 2>&1; echo "$datasetDir")
    [[ -d $d ]] || { echo "No datasetDir for $c ($d)"; exit 1; }
    dirs+=("${c%_short}=$d")
done
echo "fastANI: $($FASTANI --version 2>&1 | head -1)"
echo "pairs: $PAIRS"
echo "skani: $SKANI_OUT"

$PY - "$FASTANI" "$THREADS" "$PAIRS" "$SKANI_OUT" "$OUT" "${dirs[@]}" <<'PYEOF'
import csv, os, subprocess, sys, tempfile
fastani, threads, pairs_tsv, skani_tsv, out = sys.argv[1:6]
dirs = dict(a.split("=", 1) for a in sys.argv[6:])
KEYS = ["Cohort", "PatientID", "ARO", "gtdb_species", "SampleID_before", "MAG_before", "SampleID_after", "MAG_after"]
FLAGSHIP = {"D20248", "H04", "H75", "H61"}

def cohort(row, sid):
    if row["Cohort"] == "Heston":
        return "Heston"
    return "UNC" if sid.startswith("BMT") else "Duke"

def mag_path(row, w):
    sid = row[f"SampleID_{w}"]
    return os.path.join(dirs[cohort(row, sid)], "3.3_reassemble_bins", sid,
                        "reassembled_bins", row[f"MAG_version_{w}"] + ".fa")

rows = [r for r in csv.DictReader(open(pairs_tsv, newline=""), delimiter="\t") if r.get("tier2_accessory") == "TRUE"]
skani = {tuple(r[k] for k in KEYS): r for r in csv.DictReader(open(skani_tsv, newline=""), delimiter="\t")}

missing = sorted({p for r in rows for p in (mag_path(r, "before"), mag_path(r, "after")) if not os.path.exists(p)})
if missing:
    sys.exit(f"{len(missing)} MAG fasta(s) missing, e.g. {missing[:3]}")

def run(q, r, tmp):
    o = os.path.join(tmp, "o.txt")
    p = subprocess.run([fastani, "-q", q, "-r", r, "-t", threads, "-o", o], capture_output=True, text=True)
    if p.returncode != 0:
        return ("NA", "NA", "NA", "fastani_error")
    line = open(o).read().split()
    if not line:
        return ("NA", "NA", "NA", "below_fastani_cutoff")
    return (line[2], line[3], line[4], "")

# flagship patients first, so they appear at the top of the log
rows.sort(key=lambda r: r["PatientID"] not in FLAGSHIP)
res = {}
with tempfile.TemporaryDirectory() as tmp:
    for r in rows:
        b, a = mag_path(r, "before"), mag_path(r, "after")
        if (b, a) in res:
            continue
        fwd = run(a, b, tmp)   # query after, ref before (same direction as skani)
        rev = run(b, a, tmp)
        anis = [float(x[0]) for x in (fwd, rev) if x[0] != "NA"]
        mean = f"{sum(anis) / len(anis):.4f}" if anis else "NA"
        note = ";".join(sorted({x[3] for x in (fwd, rev) if x[3]}))
        res[(b, a)] = (fwd, rev, mean, note)
        print(f"  {len(res):>4}  {r['PatientID']:<7} {r['gtdb_species']:<35} fastANI {mean:>8}  "
              f"(fwd {fwd[0]}, rev {rev[0]})  {note}", flush=True)

def call(ani):
    if ani in ("", "NA", None):
        return "no ANI"
    a = float(ani)
    return "same strain" if a >= 99.9 else "unresolved" if a >= 99 else "different strain"

with open(out, "w", newline="") as fo:
    w = csv.writer(fo, delimiter="\t", lineterminator="\n")
    w.writerow(KEYS + ["Gene", "skani_ANI", "skani_call", "fastANI_after_vs_before", "frag_matched_fwd", "frag_total_fwd",
                       "fastANI_before_vs_after", "frag_matched_rev", "frag_total_rev", "fastANI_mean", "fastANI_call", "fastani_note"])
    tab = {}
    flag_lines = []
    for r in rows:
        fwd, rev, mean, note = res[(mag_path(r, "before"), mag_path(r, "after"))]
        s = skani.get(tuple(r[k] for k in KEYS), {}).get("ANI", "NA")
        sc, fc = call(s), call(mean)
        tab[(sc, fc)] = tab.get((sc, fc), 0) + 1
        w.writerow([r[k] for k in KEYS] + [r["Gene"], s, sc, fwd[0], fwd[1], fwd[2], rev[0], rev[1], rev[2], mean, fc, note])
        if r["PatientID"] in FLAGSHIP:
            flag_lines.append(f"  {r['PatientID']:<7} {r['Gene']:<28} {r['gtdb_species']:<30} skani {s:>8} ({sc})  fastANI {mean:>8} ({fc})")

print(f"\n{out}: {len(rows)} accessory events, {len(res)} unique MAG pairs")
print("\nFlagship events:")
print("\n".join(flag_lines))
print("\nskani call (rows) vs fastANI call (columns), accessory events:")
cats = ["same strain", "unresolved", "different strain", "no ANI"]
print(f"  {'':<18}" + "".join(f"{c:>18}" for c in cats))
for sc in cats:
    print(f"  {sc:<18}" + "".join(f"{tab.get((sc, fc), 0):>18}" for fc in cats))
PYEOF
echo "done"
