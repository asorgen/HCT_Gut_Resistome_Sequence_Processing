#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=humann_summary
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=200GB
#SBATCH --time=04:00:00

# Aggregate per-sample HUMAnN output (5.8_humann/<ID>/<ID>_{genefamilies,pathabundance}.tsv,
# 214 Duke samples) into cohort-wide tables, for the gene-family DASC correlation
# companion analysis ([[Build taxonomy and gene-family companion analyses for
# SNP paper]]). Run once, after every 5.8_humann sample has completed
# (check: `ls .../5.8_humann/COMPLETE | wc -l` should be 214).
#
# Pipeline: humann_join_tables -> humann_regroup_table (UniRef90 -> KO, using
#   the bundled map_ko_uniref90.txt.gz) -> humann_renorm_table (CPM) ->
#   humann_split_stratified_table (separate per-species-stratified rows from
#   the whole-community "unstratified" total each KO/pathway sums to - the
#   companion analysis uses UNSTRATIFIED tables, matching the taxonomy-level
#   analysis's whole-community framing).
#
# IMPORTANT (2026-07-24 incident): the joined stratified genefamilies table is
# huge (214 samples x every UniRef90 hit incl. per-species breakdown - came
# out to 4.3GB / 7,053,818 rows) and a first attempt OOM'd at 64GB trying to
# CPM-normalize it directly. Fix: regroup UniRef90 -> KO FIRST (collapses to
# ~12k rows per earlier single-sample testing, >99% size reduction) and only
# CPM-normalize the much smaller regrouped table - never renormalize the raw
# genefamilies table at full UniRef90 resolution. Bumped to 200GB regardless,
# since the regroup step itself still has to load the full 4.3GB join once.
# Each step checks for its own output first so a rerun resumes rather than
# redoing the (slow) joins.
#
# Usage: sbatch humann_summary.sh

set -euo pipefail
module load humann/3.8

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../configs/private.config"

ROOT=${HPC_PROJECTS}/HCT_Gut_Resistome_Study
HUMANN_DIR=${ROOT}/HCT_Gut_Resistome_Data/unprocessed/Duke/Duke_short/5.8_humann
OUT_DIR=${ROOT}/HCT_Gut_Resistome_Data/processed/Duke/Duke_short_tables
KO_MAP=/projects/datasets/humann/3.8/utility_mapping/map_ko_uniref90.txt.gz
mkdir -p "$OUT_DIR"

N_COMPLETE=$(ls "${HUMANN_DIR}/COMPLETE" | wc -l)
echo "HUMAnN samples complete: ${N_COMPLETE}"

GENEFAM=${OUT_DIR}/Duke_short_humann_genefamilies.tsv
PATHAB=${OUT_DIR}/Duke_short_humann_pathabundance.tsv
KO_RAW=${OUT_DIR}/Duke_short_humann_ko.tsv
KO_CPM=${OUT_DIR}/Duke_short_humann_ko_cpm.tsv
PATHAB_CPM=${OUT_DIR}/Duke_short_humann_pathabundance_cpm.tsv

if [ -s "$GENEFAM" ]; then
    echo "=== Gene-family join already found, skipping ==="
else
    echo "=== Joining gene-family tables ==="
    humann_join_tables --input "$HUMANN_DIR" --output "$GENEFAM" \
        --file_name genefamilies --search-subdirectories
fi

if [ -s "$PATHAB" ]; then
    echo "=== Pathway-abundance join already found, skipping ==="
else
    echo "=== Joining pathway-abundance tables ==="
    humann_join_tables --input "$HUMANN_DIR" --output "$PATHAB" \
        --file_name pathabundance --search-subdirectories
fi

if [ -s "$KO_RAW" ]; then
    echo "=== KO regroup already found, skipping ==="
else
    echo "=== Regrouping UniRef90 gene families -> KO (on RAW counts, before CPM - see header) ==="
    humann_regroup_table --input "$GENEFAM" \
        --output "$KO_RAW" --groups uniref90_ko --custom "$KO_MAP"
fi

if [ -s "$KO_CPM" ]; then
    echo "=== KO CPM normalization already found, skipping ==="
else
    echo "=== Normalizing KO table to CPM (small now - regrouped, not raw UniRef90) ==="
    humann_renorm_table --input "$KO_RAW" --output "$KO_CPM" --units cpm
fi

if [ -s "$PATHAB_CPM" ]; then
    echo "=== Pathway CPM normalization already found, skipping ==="
else
    echo "=== Normalizing pathway-abundance table to CPM ==="
    humann_renorm_table --input "$PATHAB" --output "$PATHAB_CPM" --units cpm
fi

echo "=== Splitting stratified/unstratified (KO, pathway) ==="
humann_split_stratified_table --input "$KO_CPM" --output "${OUT_DIR}"
humann_split_stratified_table --input "$PATHAB_CPM" --output "${OUT_DIR}"

echo "Done. Unstratified (whole-community) tables for the DASC correlation step:"
echo "  ${OUT_DIR}/Duke_short_humann_ko_cpm_unstratified.tsv"
echo "  ${OUT_DIR}/Duke_short_humann_pathabundance_cpm_unstratified.tsv"
