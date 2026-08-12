#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=uamr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32GB
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/asorgen/hct_compare/logs/uamr_%A_%a.log

# ------------------------------------------------------------------------------
# Uniform ARG calling across the short / long / hybrid assembly tracks.
#
# Why this exists rather than reusing each pipeline's own AMR output:
# the Duke_short bash pipeline and the Nanopore Snakemake pipeline do NOT call ARGs
# the same way, so their existing outputs are not directly comparable --
#   * short_scripts/5.5_AA_amr_assembly.sh runs `prodigal` with no `-p meta`;
#     the nanopore pipeline runs `prodigal -p meta` (metagenome ORF calling).
#   * short_scripts/5.5_AA_amr_assembly.sh runs `rgi load --card_json` WITHOUT
#     `--card_annotation ... --local`, the pattern that can silently drop every hit
#     when the packaged protein reference is stale relative to card.json. The short
#     calls on disk predate that fix.
# Running one identical stack over all three assemblies means the short/long/hybrid
# differences reflect the assemblies, not the tooling.
#
# All three inputs are contig-length filtered at >=1500 bp before ORF calling
# (short and long are already filtered by their pipelines; hybrid is filtered here).
#
# Usage:  sbatch --array=0-4 uniform_amr_compare.sbatch <short|long|hybrid>
#
# @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics
# ------------------------------------------------------------------------------

set -euo pipefail

TRACK=${1:?"usage: sbatch --array=0-4 uniform_amr_compare.sbatch <short|long|hybrid>"}

SAMPLES=(D20248PRE D20248D1 D20248D21 D21309PRE D21309D10)
ID=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

STUDY=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Study
CONDA_ROOT=${STUDY}/Nanopore_Resistome_Pipeline/.snakemake/conda
AMR_ENV=${CONDA_ROOT}/7aee80a0343ab5e98ae28a6ff5c374ec_
ASM_ENV=${CONDA_ROOT}/c3d6b3cd990503a915fa74b3f0d80cc2_
CARD_DB=/scratch/asorgen/RGI_databases/rgi_4.0.1

WORK=/scratch/asorgen/hct_compare
OUTDIR=${WORK}/amr/${TRACK}/${ID}
mkdir -p "${OUTDIR}" "${WORK}/logs"

# ---- resolve the assembly for this track -------------------------------------
case "${TRACK}" in
    short)
        ASSEMBLY=${STUDY}/HCT_Gut_Resistome_Data/unprocessed/Duke/Duke_short/1.2_evaluation/${ID}_final_assembly.fasta
        ;;
    long)
        ASSEMBLY=${STUDY}/Nanopore_Resistome_Pipeline/results/${ID}/evaluate/${ID}_final_assembly.fasta
        ;;
    hybrid)
        RAW=/scratch/asorgen/hct_hybrid/opera_ms/${ID}/assembly.fasta
        [[ -s "${RAW}" ]] || { echo "[ERROR] hybrid assembly missing: ${RAW}"; exit 1; }
        ASSEMBLY=${OUTDIR}/${ID}_final_assembly.fasta
        echo "[ filtering hybrid contigs to >=1500 bp ]"
        "${ASM_ENV}/bin/reformat.sh" in="${RAW}" out="${ASSEMBLY}" minlength=1500 overwrite=true
        ;;
    *)
        echo "[ERROR] unknown track: ${TRACK} (expected short|long|hybrid)"; exit 1 ;;
esac

[[ -s "${ASSEMBLY}" ]] || { echo "[ERROR] assembly not found: ${ASSEMBLY}"; exit 1; }

echo "===== ${TRACK} / ${ID} ====="
echo "[ host     ] $(hostname)"
echo "[ assembly ] ${ASSEMBLY} ($(du -h "${ASSEMBLY}" | cut -f1), $(grep -c '^>' "${ASSEMBLY}") contigs)"
date
SECONDS=0

export PATH="${AMR_ENV}/bin:${PATH}"

# ---- 1. ORF prediction (identical settings across tracks) --------------------
FAA=${OUTDIR}/${ID}_genes.faa
if [[ ! -s "${FAA}" ]]; then
    echo "[ prodigal -p meta ]"
    "${AMR_ENV}/bin/prodigal" -i "${ASSEMBLY}" -a "${FAA}" -p meta -q
fi
echo "[ ORFs ] $(grep -c '^>' "${FAA}")"

# ---- 2. AMRFinder+ (protein mode) --------------------------------------------
if [[ ! -s "${OUTDIR}/${ID}.amrfinder.tsv" ]]; then
    echo "[ amrfinder -p ]"
    "${AMR_ENV}/bin/amrfinder" --threads 8 -p "${FAA}" -o "${OUTDIR}/${ID}.amrfinder.tsv"
fi

# ---- 3. RGI (protein mode, --local against a matching CARD annotation) -------
# RGI chokes on prodigal's trailing '*' stop marker; strip it first.
# `cd` into the per-sample dir so parallel array tasks don't clobber a shared localDB/.
CLEAN_FAA=${OUTDIR}/${ID}_genes_rgi.faa
if [[ ! -s "${OUTDIR}/${ID}.rgi.txt" ]]; then
    tr -d '*' < "${FAA}" > "${CLEAN_FAA}"
    clean_abs=$(realpath "${CLEAN_FAA}")
    cd "${OUTDIR}"
    echo "[ rgi load --local ]"
    "${AMR_ENV}/bin/rgi" load --card_json "${CARD_DB}/card.json" \
        --card_annotation "${CARD_DB}/card_database_v4.0.1.fasta" --local
    echo "[ rgi main -t protein ]"
    "${AMR_ENV}/bin/rgi" main -i "$clean_abs" -t protein -o "${ID}.rgi" --local --clean -n 8
fi

# ---- report ------------------------------------------------------------------
af_n=$(($(wc -l < "${OUTDIR}/${ID}.amrfinder.tsv") - 1))
rgi_n=$(($(wc -l < "${OUTDIR}/${ID}.rgi.txt") - 1))
echo "[ AMRFinder+ calls ] ${af_n}"
echo "[ RGI calls        ] ${rgi_n}"
if [[ "${rgi_n}" -eq 0 ]]; then
    echo "[ WARNING ] RGI returned 0 calls -- check for a stale CARD protein reference."
fi

touch "${OUTDIR}/COMPLETE"
echo "[ elapsed ] $((SECONDS / 60))m"
echo "===== ${TRACK} / ${ID} DONE ====="
