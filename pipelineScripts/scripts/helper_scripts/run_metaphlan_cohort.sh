#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96GB
#SBATCH --time=02:00:00
# 96GB: MaxRSS sat at 61-64 GB against the old 64 GB request (Duke run + 2026-11-10 smoke test)

# Standalone MetaPhlAn4 (module 2.3) submission for one short-read cohort.
#
# Runs the unchanged 2.3_metaphlan4.sh as a throttled SLURM array, outside
# run_short_pipeline.sh, so that enabling MetaPhlAn does not also resubmit
# whatever other modules the cohort config has switched on (e.g. Heston's
# PHASE 4 binning). Outputs and COMPLETE flags land in the same
# ${datasetDir}/2.3_metaphlan4 the wrapper would use, so a later wrapper run
# with run_m4=true sees them as done.
#
# Usage (from anywhere, on the login node — it only builds a list and submits):
#   bash run_metaphlan_cohort.sh <config> [max_concurrent]
#     config          e.g. Duke_short, UNC_short, Heston_short
#     max_concurrent  array throttle (default 40)
#   MAX_SAMPLES=1 bash run_metaphlan_cohort.sh Duke_short   # smoke test
#
# Samples = those with 0.4 decontaminated reads, minus read-depth EXCLUDED,
# minus samples already COMPLETE for 2.3.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # ---------------- Submission mode ----------------
    CONFIG_NAME=$1
    MAX_CONC=${2:-40}
    CONFIG_DIR="$(cd "${SCRIPT_DIR}/../../configs" && pwd)"
    CONFIG_FILE=${CONFIG_DIR}/${CONFIG_NAME}-read.config
    [[ -f "$CONFIG_FILE" ]] || { echo "No such config: $CONFIG_FILE"; exit 1; }

    source "$CONFIG_FILE"
    [[ -n "$datasetDir" && -n "$metaphlanDir" && -n "$clean_readDir" ]] || { echo "Config missing datasetDir/metaphlanDir/clean_readDir"; exit 1; }

    outDir=${metaphlanDir}
    mkdir -p ${outDir}/logs ${outDir}/COMPLETE
    idList=${outDir}/ids_$(date +%Y%m%d_%H%M%S).txt

    excl=${datasetDir}/0.5_read_depth_filter/EXCLUDED
    for r1 in ${clean_readDir}/*_1.fastq.gz; do
        ID=$(basename "$r1" _1.fastq.gz)
        [[ -e ${excl}/${ID} ]] && continue
        [[ -e ${outDir}/COMPLETE/${ID} ]] && continue
        [[ -s ${clean_readDir}/${ID}_2.fastq.gz ]] || { echo "WARN: no R2 for ${ID}, skipping"; continue; }
        echo "$ID"
    done > "$idList"

    # Optional cap (smoke tests): keep only the first N samples
    if [[ -n "${MAX_SAMPLES:-}" ]]; then
        head -n "$MAX_SAMPLES" "$idList" > "${idList}.tmp" && mv "${idList}.tmp" "$idList"
    fi

    N=$(wc -l < "$idList")
    echo "${CONFIG_NAME}: ${N} samples to run -> ${idList}"
    [[ $N -eq 0 ]] && { rm -f "$idList"; exit 0; }

    sbatch --array=1-${N}%${MAX_CONC} \
        --chdir=${datasetDir} \
        --job-name=2.3_m4_${dataset} \
        -o ${outDir}/logs/%x.%A_%a.log \
        --export=ALL,M4_CONFIG=${CONFIG_FILE},M4_IDLIST=${idList} \
        "${SCRIPT_DIR}/run_metaphlan_cohort.sh"
    exit $?
fi

# ---------------- Array-task mode ----------------
# Mirror the environment run_short_pipeline.sh + module_setup give a module script.
source "$M4_CONFIG"
export pipelineConfig=$M4_CONFIG
export config_file=$(which config-metawrap 2>/dev/null || echo /dev/null)
export bashrc bash_profile clean_readDir metaphlanDir METAPHLAN4_DB METAPHLAN4_INDEX
export ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$M4_IDLIST")
[[ -n "$ID" ]] || { echo "No ID at line ${SLURM_ARRAY_TASK_ID} of ${M4_IDLIST}"; exit 1; }
export moduleDir=2.3_metaphlan4   # relative to --chdir=${datasetDir}, as module_setup does

cd "$datasetDir" || exit 1
bash ${moduleScriptPath}/2.3_metaphlan4.sh
rc=$?

# --mapout is ~180 MB/sample uncompressed (~200 GB for 1,110 samples on 94-97% full
# filesystems). Keep it (StrainPhlAn input) but bzip2 it; MetaPhlAn reads .bz2 natively.
map=${moduleDir}/${ID}_mapping.txt
if [[ $rc -eq 0 && -s $map ]]; then bzip2 -f "$map"; fi
exit $rc
