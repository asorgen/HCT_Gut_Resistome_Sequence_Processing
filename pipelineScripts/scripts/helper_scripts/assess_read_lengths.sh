#!/bin/bash
#SBATCH --job-name=assess_read_lengths
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --time=8:00:00
#SBATCH --output=/projects/afodor_research/asorgen/datasets/PRJNA890666/Heston_short/read_length_assessment/assess_%j.out
#SBATCH --mail-type=END,FAIL

# Measure the ACTUAL read-length distribution of post-QC reads.
#
#   sbatch assess_read_lengths.sh <dataset_dir> [module_subdir]
#
# Why post-QC rather than SRA metadata: PRJNA890666 was submitted already
# quality-trimmed with variable read lengths (FASTQ headers still claim
# length=151 while reads measure ~130), and this pipeline's own trimming shifts
# them further. 2.1_kraken2.sh picks the Bracken kmer_distrib from the length of
# the reads it is actually given, which are these files - so this is the
# distribution that determines the 100mer/150mer split.
#
# Emits per-file stats (min/avg/max/Q1/Q2/Q3/N50) for every sample.

set -o pipefail

DATASET_DIR=${1:-/projects/afodor_research/asorgen/datasets/PRJNA890666/Heston_short}
MODULE_SUBDIR=${2:-0.4_host_decontamination}

READ_DIR=${DATASET_DIR}/${MODULE_SUBDIR}
OUT_DIR=${DATASET_DIR}/read_length_assessment
mkdir -p "$OUT_DIR"

echo "[$(date)] read-length assessment"
echo "  reads : ${READ_DIR}"
echo "  out   : ${OUT_DIR}"
echo "  host  : $(hostname), ${SLURM_CPUS_PER_TASK} cpus"

n_files=$(ls "${READ_DIR}"/*.fastq.gz 2>/dev/null | wc -l)
if [[ $n_files -eq 0 ]]; then
    echo "[$(date)] FAILED - no FASTQs in ${READ_DIR}"
    exit 1
fi
echo "  files : ${n_files}"

module load seqkit/2.9.0

SECONDS=0
# -a gives the full stat set (Q1/Q2/Q3/N50), -T tab-delimited for downstream parsing.
seqkit stats -a -T -j "${SLURM_CPUS_PER_TASK}" "${READ_DIR}"/*.fastq.gz \
    > "${OUT_DIR}/seqkit_stats_raw.tsv"
STATUS=$?

if [[ $STATUS -ne 0 || ! -s "${OUT_DIR}/seqkit_stats_raw.tsv" ]]; then
    echo "[$(date)] FAILED - seqkit exited ${STATUS}"
    exit 1
fi

echo "[$(date)] seqkit done in ${SECONDS}s ($(wc -l < "${OUT_DIR}/seqkit_stats_raw.tsv") rows)"
echo "[$(date)] SUCCESS - ${OUT_DIR}/seqkit_stats_raw.tsv"
