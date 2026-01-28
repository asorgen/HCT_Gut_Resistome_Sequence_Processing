#!/bin/bash

# Shell script to automate MetaPhlAn2 
# Author = Alicia Sorgen
# Date = 2025 Nov 3


# Set up

input=""
output=""
module_dir=""
# inputList=${input}/BariatricSurgery_FileList.txt

# Help message
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -i, --input     Path to merged MetaPhlAn4 input file (required)
  -o, --output    Path to output directory and pipeline name (required)
  -p, --profiles  Path to directory containing MetaPhlAn4 profile files (required)
  -h, --help      Show this help message and exit

Example:
  $(basename "$0") -i data/input.txt -o results/output.txt
EOF
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            input="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
            shift 2
            ;;
        -p|--profiles)
            module_dir="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown argument: $1"
            usage
            ;;
    esac
done

# Validate required args
if [[ -z "$input" || -z "$output" ]]; then
    echo "Error: --input, --output, and --profiles are required."
    usage
fi

# Main logic
echo "Input: $input" # ${out}/${pipeline}_metaphlan4_counts.tsv
echo "Output: $output" # ${out}/${pipeline}
echo "Profile directory: $module_dir"


funcScript=${ps_path}/functions.R

# export input
# export output
# export funcScript
# export ps_path

### Load modules
module load R/4.3.3 

Rscript \
${ps_path}/BuildTaxaTable_MetaPhlAn4.R \
${input} \
${output} \
${funcScript} ${module_dir}

# Rscript \
# ${ps_path}/NormalizeCounts_MetaPhlAn2.R \
# ${input} \
# ${output} \
# ${funcScript}

