#!/bin/bash

# annotate_bins_summary.sh annotate_bins bin_funct_annotations SHORT_bin_classification.tsv

pythonScript="scripts/post_scripts/annotate_bins_summary.py"

# Set the parent directory
PARENT_DIR=$1
ANNOTATION_DIR=$2
MASTER_FILE=${PARENT_DIR}/$3

# Loop through each directory within the parent directory
for dir in "$PARENT_DIR"/*/; do
  
  # Extract only the directory name (basename)
  dir_name=$(basename "$dir")

  if [[ -d $dir/${ANNOTATION_DIR} ]]; then
    for file in `ls ${dir}/${ANNOTATION_DIR}`; do
      file_name=${dir}${ANNOTATION_DIR}/${file}
      # echo $file_name
      python $pythonScript "$file_name" "$dir_name" "$MASTER_FILE"
    done
  fi

  
done


