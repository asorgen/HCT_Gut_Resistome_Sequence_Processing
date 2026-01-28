#!/bin/bash

# gene_profiling_summary.sh asm_gene_profiling gene_annotations.tsv SHORT_gene_annotations.tsv

pythonScript="${ps_path}/gene_profiling_summary.py"

# Set the parent directory
PARENT_DIR=$1
INPUT_FILE=$2
MASTER_FILE=${PARENT_DIR}/$3

# Loop through each directory within the parent directory
for dir in "$PARENT_DIR"/*; do

  # Extract only the directory name (basename)
  dir_name=$(basename "$dir") 

  file=${dir}/${dir_name}_${INPUT_FILE} 

  if [[ -f $file ]]; then
    python $pythonScript "$file" "$dir_name" "$MASTER_FILE"
  fi

done