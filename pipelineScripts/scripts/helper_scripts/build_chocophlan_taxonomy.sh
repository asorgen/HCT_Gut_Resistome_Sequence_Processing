#!/bin/bash
#SBATCH --job-name=build_chocophlan_taxonomy
#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=2:00:00
#SBATCH --output=/scratch/asorgen/waafle_db/build_chocophlan_taxonomy_%j.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aliciasorgen0@gmail.com

# One-time setup: build the WAAFLE taxonomy TSV from ChocoPhlAn FASTA headers.
#
# WAAFLE's waafle_orgscorer requires a taxonomy file mapping each clade to its
# parent (2-column TSV: clade<TAB>parent). Lineages are embedded in the
# ChocoPhlAn sseqid as field [1] of the '|'-delimited header (e.g.,
# "k__Bacteria.p__Firmicutes....s__Enterococcus_faecium"). This script reads
# all .ffn.gz headers, extracts unique lineage strings, and writes all
# clade->parent pairs up to the WAAFLE root ("r__Root").
#
# Output: /scratch/asorgen/waafle_db/chocophlan.taxonomy
# Created: 2026-05-26

CHOCOPHLAN_SRC=/projects/datasets/humann/3.8/chocophlan
TAXONOMY_OUT=/scratch/asorgen/waafle_db/chocophlan.taxonomy

echo "[$(date)] Building WAAFLE taxonomy from ChocoPhlAn headers"
echo "Source: $CHOCOPHLAN_SRC"
echo "Output: $TAXONOMY_OUT"

module load anaconda3/2023.09
source /users/asorgen/miniforge3/etc/profile.d/conda.sh
conda activate waafle-env

python3 - <<'PYEOF'
import gzip, os, sys

chocophlan_src = "/projects/datasets/humann/3.8/chocophlan"
taxonomy_out   = "/scratch/asorgen/waafle_db/chocophlan.taxonomy"
root_clade     = "r__Root"

pairs = set()   # (clade, parent) pairs

ffn_files = [f for f in os.listdir(chocophlan_src) if f.endswith(".ffn.gz")]
print(f"[build] Reading {len(ffn_files)} .ffn.gz files ...", flush=True)

seen_lineages = set()
for i, fname in enumerate(ffn_files):
    if i % 1000 == 0:
        print(f"  {i}/{len(ffn_files)}", flush=True)
    path = os.path.join(chocophlan_src, fname)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip()
            fields = header.split("|")
            if len(fields) < 2:
                continue
            lineage = fields[1]   # e.g. k__Bacteria.p__Firmicutes....s__Enterococcus_faecium
            if lineage in seen_lineages:
                continue
            seen_lineages.add(lineage)
            # Build clade->parent pairs by progressively stripping last level
            levels = lineage.split(".")
            for depth in range(len(levels), 0, -1):
                clade  = ".".join(levels[:depth])
                parent = ".".join(levels[:depth - 1]) if depth > 1 else root_clade
                pairs.add((clade, parent))

# Add kingdom-level -> root pairs already covered above, but ensure root itself present
print(f"[build] Unique lineages: {len(seen_lineages)}", flush=True)
print(f"[build] Taxonomy pairs:  {len(pairs)}", flush=True)

with open(taxonomy_out, "w") as fh:
    for clade, parent in sorted(pairs):
        fh.write(f"{clade}\t{parent}\n")

print(f"[build] Written to {taxonomy_out}", flush=True)
PYEOF

STATUS=$?
if [[ $STATUS -eq 0 ]]; then
    echo "[$(date)] SUCCESS"
    wc -l "$TAXONOMY_OUT"
else
    echo "[$(date)] FAILED (exit $STATUS)"
    exit 1
fi
