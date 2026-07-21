#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --job-name=humann_summary
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64GB
#SBATCH --time=04:00:00

# Aggregate per-sample HUMAnN output (5.8_humann/<ID>/<ID>_{genefamilies,pathabundance}.tsv,
# 214 Duke samples) into cohort-wide tables, for the gene-family DASC correlation
# companion analysis ([[Build taxonomy and gene-family companion analyses for
# SNP paper]]). Run once, after every 5.8_humann sample has completed
# (check: `ls .../5.8_humann/COMPLETE | wc -l` should be 214).
#
# Pipeline: humann_join_tables -> humann_renorm_table (CPM) ->
#   humann_regroup_table (UniRef90 -> KO, using the bundled map_ko_uniref90.txt.gz) ->
#   humann_split_stratified_table (separate per-species-stratified rows from
#   the whole-community "unstratified" total each KO/pathway sums to - the
#   companion analysis uses UNSTRATIFIED tables, matching the taxonomy-level
#   analysis's whole-community framing).
#
# Usage: sbatch humann_summary.sh

set -euo pipefail
module load humann/3.8

ROOT=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Study
HUMANN_DIR=${ROOT}/HCT_Gut_Resistome_Data/unprocessed/Duke/Duke_short/5.8_humann
OUT_DIR=${ROOT}/HCT_Gut_Resistome_Data/processed/Duke/Duke_short_tables
KO_MAP=/projects/datasets/humann/3.8/utility_mapping/map_ko_uniref90.txt.gz
mkdir -p "$OUT_DIR"

N_COMPLETE=$(ls "${HUMANN_DIR}/COMPLETE" | wc -l)
echo "HUMAnN samples complete: ${N_COMPLETE}"

echo "=== Joining gene-family tables ==="
humann_join_tables --input "$HUMANN_DIR" --output "${OUT_DIR}/Duke_short_humann_genefamilies.tsv" \
    --file_name genefamilies --search-subdirectories

echo "=== Joining pathway-abundance tables ==="
humann_join_tables --input "$HUMANN_DIR" --output "${OUT_DIR}/Duke_short_humann_pathabundance.tsv" \
    --file_name pathabundance --search-subdirectories

echo "=== Normalizing to CPM ==="
humann_renorm_table --input "${OUT_DIR}/Duke_short_humann_genefamilies.tsv" \
    --output "${OUT_DIR}/Duke_short_humann_genefamilies_cpm.tsv" --units cpm
humann_renorm_table --input "${OUT_DIR}/Duke_short_humann_pathabundance.tsv" \
    --output "${OUT_DIR}/Duke_short_humann_pathabundance_cpm.tsv" --units cpm

echo "=== Regrouping UniRef90 gene families -> KO ==="
humann_regroup_table --input "${OUT_DIR}/Duke_short_humann_genefamilies_cpm.tsv" \
    --output "${OUT_DIR}/Duke_short_humann_ko_cpm.tsv" --groups uniref90_ko --custom "$KO_MAP"

echo "=== Splitting stratified/unstratified (KO, pathway) ==="
humann_split_stratified_table --input "${OUT_DIR}/Duke_short_humann_ko_cpm.tsv" --output "${OUT_DIR}"
humann_split_stratified_table --input "${OUT_DIR}/Duke_short_humann_pathabundance_cpm.tsv" --output "${OUT_DIR}"

echo "Done. Unstratified (whole-community) tables for the DASC correlation step:"
echo "  ${OUT_DIR}/Duke_short_humann_ko_cpm_unstratified.tsv"
echo "  ${OUT_DIR}/Duke_short_humann_pathabundance_cpm_unstratified.tsv"
