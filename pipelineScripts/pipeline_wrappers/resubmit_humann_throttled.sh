#!/bin/bash
# Throttled resubmission for 5.8_humann.sh failures (2026-07-21 disk-space
# incident - see that script's header comment). Submits one sample at a time
# via test_single_sample.sh, capping concurrent HUMAnN jobs in queue/running
# so scratch staging usage (HPC_SCRATCH/humann_staging/, ~tens of GB per
# sample while running) stays well within headroom, rather than firing all
# failures at once like the original full-cohort submission did.
#
# Usage: bash resubmit_humann_throttled.sh <failed_ids_file> [max_concurrent]
#   failed_ids_file: one sample ID per line (e.g. from `comm -23` against
#                     the COMPLETE dir - see incident notes)
#   max_concurrent:  default 20

set -eo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

IDS_FILE="${1:?Usage: resubmit_humann_throttled.sh <failed_ids_file> [max_concurrent]}"
MAX_CONCURRENT="${2:-20}"

n_total=$(wc -l < "$IDS_FILE")
n_submitted=0

echo "=== Throttled HUMAnN resubmission: ${n_total} samples, max ${MAX_CONCURRENT} concurrent ==="

while IFS= read -r id; do
    [ -z "$id" ] && continue
    while true; do
        # grep -c exits 1 (not just "0") when nothing matches, which would
        # trip `set -e` here - `|| true` makes the zero-match case explicit
        # instead of silently killing the whole resubmission loop.
        n_queued=$(squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -c '^5\.8_' || true)
        if [ "$n_queued" -lt "$MAX_CONCURRENT" ]; then
            break
        fi
        sleep 60
    done
    bash "${SCRIPT_DIR}/test_single_sample.sh" Duke short 5.8_humann.sh "$id" humann_opts
    n_submitted=$((n_submitted + 1))
    echo "[$n_submitted/$n_total] submitted $id"
done < "$IDS_FILE"

echo "=== All ${n_submitted} resubmissions queued ==="
