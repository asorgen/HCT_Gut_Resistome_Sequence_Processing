#!/bin/bash
# Throttled driver for run_short_pipeline.sh.
#
# Re-invokes the wrapper with MAX_SAMPLES_PER_RUN, waiting for the queue to drain
# below a threshold between invocations. Needed because the wrapper's cap limits
# submissions PER INVOCATION, not total in-flight — calling it again while the
# previous samples are still queued would stack another batch on top.
#
# Why throttle at all:
#   1. Orion's QOS caps submitted jobs at 2048, and the wrapper submits one job
#      per module per sample (754 samples x 5 modules = 3770 for Heston QC).
#   2. Peak disk is driven by in-flight samples, not cohort size — 0.1_pre_qc
#      stages uncompressed FASTQ (~8.5 GB/sample for Heston) that is only
#      reclaimed by the cleanup module at the end of each sample's chain.
#
#   Usage: bash run_short_pipeline_throttled.sh <config> [batch_size] [drain_to]
#     config      config name (e.g. Heston_short) or path to a *-read.config
#     batch_size  samples to submit per invocation (default 25)
#     drain_to    resubmit once in-flight jobs for this dataset fall to or below
#                 this number (default 0, i.e. fully drained)
#
#   Honours SAMPLE_LIST_OVERRIDE if exported.
#
# Safe to interrupt and restart: progress lives in the per-module COMPLETE flags,
# so a restart picks up wherever it left off.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

CONFIG_ARG=${1:?Usage: run_short_pipeline_throttled.sh <config> [batch_size] [drain_to]}
BATCH_SIZE=${2:-25}
DRAIN_TO=${3:-0}
POLL=60

cd "$PIPELINE_ROOT" || exit 1

if [[ -f "$CONFIG_ARG" ]]; then
    CONFIG_FILE=$CONFIG_ARG
else
    CONFIG_FILE=pipelineScripts/configs/${CONFIG_ARG}-read.config
fi
[[ -f "$CONFIG_FILE" ]] || { echo "No such config: ${CONFIG_FILE}"; exit 1; }

# Read dataset/sampleList/datasetDir without leaking config vars into this shell.
read -r dataset datasetDir sampleList < <(
    source "$CONFIG_FILE" >/dev/null 2>&1
    echo "$dataset" "$datasetDir" "${SAMPLE_LIST_OVERRIDE:-$sampleList}"
)
[[ -n "$dataset" ]] || { echo "Config did not define dataset"; exit 1; }

export version=${version:-$(date +"%Y.%m.%d")}
LOGDIR=${datasetDir}/LOGs
mkdir -p "$LOGDIR"
DRIVER_LOG=${LOGDIR}/${dataset}_throttled_${version}.out

total=$(grep -vc '^#' "$sampleList")

# In-flight jobs for THIS dataset only. The wrapper names jobs
# <module>_<ID>_<dataset>[...], so match on the dataset suffix to avoid counting
# other cohorts or unrelated work.
inflight() {
    squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -c "_${dataset}\(_\|$\)" || true
}

# Samples with a top-level COMPLETE flag (written by the cleanup module).
done_count() {
    ls "${datasetDir}/COMPLETE" 2>/dev/null | grep -vc '^logs$' || true
}

# Record the PID both in the log and in a pidfile, so a driver launched with
# nohup days ago can still be found and stopped without having kept the PID.
#   pgrep -af run_short_pipeline_throttled.sh
#   cat <datasetDir>/LOGs/<dataset>_throttled.pid
PIDFILE=${LOGDIR}/${dataset}_throttled.pid
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
    echo "A throttled driver for ${dataset} is already running (pid $(cat "$PIDFILE"))."
    echo "Stop it first, or delete ${PIDFILE} if it is stale."
    exit 1
fi
echo $$ > "$PIDFILE"
trap 'rm -f "$PIDFILE"' EXIT

echo "=== Throttled run: ${dataset} ===" | tee -a "$DRIVER_LOG"
echo "  pid        : $$"                 | tee -a "$DRIVER_LOG"
echo "  config     : ${CONFIG_FILE}"      | tee -a "$DRIVER_LOG"
echo "  samples    : ${total}"            | tee -a "$DRIVER_LOG"
echo "  batch size : ${BATCH_SIZE}"       | tee -a "$DRIVER_LOG"
echo "  drain to   : ${DRAIN_TO}"         | tee -a "$DRIVER_LOG"
echo "  driver log : ${DRIVER_LOG}"       | tee -a "$DRIVER_LOG"
echo "  pidfile    : ${PIDFILE}"          | tee -a "$DRIVER_LOG"

round=0
while true; do
    # Wait for the previous batch to drain before submitting more.
    while [[ $(inflight) -gt $DRAIN_TO ]]; do sleep "$POLL"; done

    before=$(done_count)
    round=$((round + 1))
    echo "[$(date '+%F %T')] round ${round}: ${before}/${total} complete, submitting up to ${BATCH_SIZE}" \
        | tee -a "$DRIVER_LOG"

    out=$(MAX_SAMPLES_PER_RUN=$BATCH_SIZE \
          bash "${SCRIPT_DIR}/run_short_pipeline.sh" "$CONFIG_ARG" 2>&1)
    echo "$out" >> "$DRIVER_LOG"

    submitted=$(echo "$out" | grep -c "^Submitted batch job" || true)
    echo "[$(date '+%F %T')] round ${round}: submitted ${submitted} jobs" | tee -a "$DRIVER_LOG"

    # Nothing submitted and nothing running means there is no work left. This is
    # the only exit condition — it covers both "all done" and "everything else is
    # permanently failing", so check the log rather than assuming success.
    if [[ $submitted -eq 0 ]]; then
        sleep 5
        if [[ $(inflight) -eq 0 ]]; then
            echo "[$(date '+%F %T')] no work submitted and queue empty - stopping" \
                | tee -a "$DRIVER_LOG"
            break
        fi
    fi
done

final=$(done_count)
echo "=== Finished: ${final}/${total} samples complete ===" | tee -a "$DRIVER_LOG"
[[ $final -lt $total ]] && \
    echo "    $((total - final)) incomplete - check module COMPLETE dirs and logs" \
        | tee -a "$DRIVER_LOG"
exit 0
