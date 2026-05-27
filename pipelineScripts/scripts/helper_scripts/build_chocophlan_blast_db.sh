#!/bin/bash
#SBATCH --job-name=build_chocophlan_db
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=4:00:00
#SBATCH --output=/scratch/asorgen/waafle_db/build_chocophlan_db_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aliciasorgen0@gmail.com

# One-time setup: concatenate ChocoPhlAn pangenome .ffn.gz files and build a
# BLAST nucleotide database for WAAFLE. Takes ~30-60 min on Orion.
#
# Output: /scratch/asorgen/waafle_db/chocophlan_blast_db.*
# Created: 2026-05-26

CHOCOPHLAN_SRC=/projects/datasets/humann/3.8/chocophlan
WAAFLE_DB_DIR=/scratch/asorgen/waafle_db
CONCAT_FASTA=${WAAFLE_DB_DIR}/chocophlan_concat.fna
BLAST_DB=${WAAFLE_DB_DIR}/chocophlan_blast_db

mkdir -p "$WAAFLE_DB_DIR"

echo "[$(date)] Starting ChocoPhlAn BLAST DB build"
echo "Source: $CHOCOPHLAN_SRC"
echo "Output: $BLAST_DB"

# Concatenate all .ffn.gz files (skip if FASTA already exists from a prior run)
if [[ -s "$CONCAT_FASTA" ]]; then
    echo "[$(date)] Concatenated FASTA already exists — skipping concat step"
else
    echo "[$(date)] Concatenating $(ls ${CHOCOPHLAN_SRC}/*.ffn.gz | wc -l) .ffn.gz files..."
    zcat ${CHOCOPHLAN_SRC}/*.ffn.gz > "$CONCAT_FASTA"
    echo "[$(date)] Concatenated FASTA: $(wc -l < $CONCAT_FASTA) lines"
fi

# Load BLAST and build database
# Note: -parse_seqids omitted — ChocoPhlAn IDs exceed the 50-char limit.
# WAAFLE uses tabular BLAST output only; sequence retrieval by ID is not needed.
module load blast/2.15.0+
echo "[$(date)] Running makeblastdb..."
makeblastdb \
    -in "$CONCAT_FASTA" \
    -dbtype nucl \
    -out "$BLAST_DB"

STATUS=$?

if [[ $STATUS -eq 0 ]]; then
    rm -f "$CONCAT_FASTA"
    echo "[$(date)] SUCCESS — BLAST database built at: $BLAST_DB"
    ls -lh ${BLAST_DB}.*
else
    echo "[$(date)] FAILED — makeblastdb exited with status $STATUS (concat FASTA preserved at $CONCAT_FASTA)"
    exit 1
fi
