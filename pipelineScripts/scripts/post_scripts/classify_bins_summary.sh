#!/bin/bash

# classify_bins_summary.sh classify_bins bin_taxonomy.tab SHORT_bin_classification.tsv

pythonScript="scripts/post_scripts/classify_bins_summary.py"

# Set the parent directory
PARENT_DIR=$1
MASTER_FILE=${PARENT_DIR}/$3

# Loop through each directory within the parent directory
for dir in "$PARENT_DIR"/*/; do
  
  # Extract only the directory name (basename)
  dir_name=$(basename "$dir")

  file=${dir}/$2
  
  # Add your commands here to process each directory
  python $pythonScript "$file" "$dir_name" "$MASTER_FILE"

done