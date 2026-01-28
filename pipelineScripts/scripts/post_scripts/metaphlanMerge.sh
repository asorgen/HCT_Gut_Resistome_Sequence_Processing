#!/bin/bash

# Shell script to merge MetaPhlAn4 data 
# Author = Alicia Sorgen
# Date = 2025 Nov 3

# Set up

module_dir=""
output=""
# inputList=${input}/BariatricSurgery_FileList.txt

# Help message
usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  -i, --input     Path to input file or directory (required)
  -o, --output    Path to output file or directory (required)
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
            module_dir="$2"
            shift 2
            ;;
        -o|--output)
            output="$2"
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
if [[ -z "$module_dir" || -z "$output" ]]; then
    echo "Error: both --input and --output are required."
    usage
fi

# Main logic
echo "Input: $module_dir" # 2.3_metaphlan4
echo "Output: $output" # ${out}/${pipeline}



if [[ ! -d ${module_dir}/bowtie2_mapping ]]; then mkdir -p ${module_dir}/bowtie2_mapping; fi
mv ${module_dir}/*_mapping.txt ${module_dir}/bowtie2_mapping

module purge
module load metaphlan/4.2.2

# if [[ -s ${output}_metaphlan4_counts.tsv ]]; then rm ${output}_metaphlan4_counts.tsv; fi
# merge_metaphlan_tables.py ${module_dir}/*_rel_abun.txt > ${output}_metaphlan4_counts.tsv

if [[ -s ${output}_metaphlan4_rel_abun.tsv ]]; then rm ${output}_metaphlan4_rel_abun.tsv; fi
merge_metaphlan_tables.py ${module_dir}/*_rel_abun_w_read_stats.txt > ${output}_metaphlan4_rel_abun.tsv

module unload metaphlan/4.2.2






