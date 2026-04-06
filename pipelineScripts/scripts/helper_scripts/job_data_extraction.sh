#!/bin/bash

log_dir=${1}logs  # change this to your log directory
echo -e "log_dir=$log_dir"
# log_dir=Duke_short/0.2_deduplication/logs

dataset="${log_dir%%/*}"
echo -e "dataset=$dataset"

module="${log_dir#*/}"
module="${module%%/*}"
echo -e "module=$module"


output_file=${dataset}/LOGs/${module}_jobs.tsv
echo -e "output_file=$output_file"

# Write header to output file
echo -e "SampleID\tJobID\tWallTime\tMemUsed\tNodes\tCores\tCPU_Util\tState\tLog" > "$output_file"


for log_file in "$log_dir"/*.log; do
    # Remove directory and file extension
    base_name=$(basename "$log_file" .log)

    # Extract JobID using parameter expansion (everything after the first dot)
    sample_id="${base_name%%.*}"
    job_id="${base_name#*.}"

    # Run seff and extract desired fields
    seff_output=$(seff "$job_id")

    state=$(echo "$seff_output" | awk -F': ' '/State:/ {print $2}' | cut -d' ' -f1)
    nodes=$(echo "$seff_output" | awk -F': ' '/Nodes:/ {print $2}')
    cores=$(echo "$seff_output" | awk -F': ' '/Cores per node:/ {print $2}')
    cpu_utilized=$(echo "$seff_output" | awk -F': ' '/CPU Utilized:/ {print $2}')
    wall_time=$(echo "$seff_output" | awk -F': ' '/Job Wall-clock time:/ {print $2}')
    mem_used=$(echo "$seff_output" | awk -F': ' '/Memory Utilized:/ {print $2}')

    # Print tab-separated output
    echo -e "$sample_id\t$job_id\t$wall_time\t$mem_used\t$nodes\t$cores\t$cpu_utilized\t$state\t${log_file}" >> "$output_file"
done

exit 0