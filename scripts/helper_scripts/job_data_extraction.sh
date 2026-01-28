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


IDs=(D13004D15
D13004D384
D13004D46
D13004D7
D13004PRE
D17010D15
D17010D2
D17010D22
D17010D31
D17010D398
D17010D64
D17010D8
D17010PRE
D18423D109
D18423D180
D18423D19
D18423D28
D18423D354
D18423D5
D18423D92
D18423PRE
D18594D5
D18594PRE
D19688D14
D19688D180
D19688D20
D19688D29
D19688D364
D19688D5
D19688D63
D19688D96
D19688PRE
D19754D110
D19754D-1
D19754D12
D19754D171
D19754D26
D19754D384
D19754D55
D19754D8
D19754D-9
D19836D0
D19836D13
D19836D189
D19836D21
D19836D365
D19836D56
D19836D92
D19836PRE
D19840D106
D19840D12
D19840D19
D19840D196
D19840D368
D19840D63
D19840D767
D19840PRE
D20117D108
D20117D12
D20117D196
D20117D20
D20117D267
D20117D29
D20117D357
D20117D4
D20117D44
D20117D537
D20117D743
D20117D76
D20117PRE
D20215D0
D20215D16
D20215D177
D20215D22
D20215D28
D20215D37
D20215D70
D20215D91
D20215PRE
D20231D1
D20231D14
D20231D188
D20231D21
D20231D29
D20231D69
D20231PRE
D20236D0
D20236D14
D20236D169
D20236D27
D20236D34
D20236D43
D20236D99)

raw_readDir=Duke_short/0.0_raw_reads
pre_qcDir=Duke_short/0.1_pre_qc
dedup_Dir=Duke_short/0.2_deduplication
trimmed_Dir=Duke_short/0.3_sequence_trim
clean_readDir=Duke_short/0.4_host_decontamination
krakenDir=Duke_short/2.1_kraken
brackenDir=Duke_short/2.2_bracken
shortbredDir=Duke_short/5.3_shortbred

for ID in "${IDs[@]}"; do
    
    # echo $ID
    Complete_tag=()

    # Complete_tag+=(${raw_readDir}/${ID}_1.fastq ${raw_readDir}/${ID}_2.fastq)
    # Complete_tag+=(${pre_qcDir}/${ID}_1_fastqc.html ${pre_qcDir}/${ID}_1_fastqc.zip ${pre_qcDir}/${ID}_2_fastqc.html ${pre_qcDir}/${ID}_2_fastqc.zip)

    Complete_tag+=(${dedup_Dir}/${ID}_deduped_R1.fastq.gz ${dedup_Dir}/${ID}_deduped_R1.fastq.gz)
    Complete_tag+=(${dedup_Dir}/${ID}_stats.log ${dedup_Dir}/QC_report/${ID}_deduped_R1_fastqc.html ${dedup_Dir}/QC_report/${ID}_deduped_R1_fastqc.zip ${dedup_Dir}/QC_report/${ID}_deduped_R2_fastqc.html ${dedup_Dir}/QC_report/${ID}_deduped_R2_fastqc.zip)

    # Complete_tag+=(${trimmed_Dir}/${ID}_trimmed_1.fastq.gz ${trimmed_Dir}/${ID}_trimmed_2.fastq.gz)
    # Complete_tag+=(${trimmed_Dir}/QC_report/${ID}_trimmed_1_fastqc.html ${trimmed_Dir}/QC_report/${ID}_trimmed_1_fastqc.zip ${trimmed_Dir}/QC_report/${ID}_trimmed_2_fastqc.html ${trimmed_Dir}/QC_report/${ID}_trimmed_2_fastqc.zip)

    # Complete_tag+=(${clean_readDir}/QC_report/${ID}_cleaned_R1_fastqc.html ${clean_readDir}/QC_report/${ID}_cleaned_R1_fastqc.zip ${clean_readDir}/QC_report/${ID}_cleaned_R2_fastqc.html ${clean_readDir}/QC_report/${ID}_cleaned_R2_fastqc.zip ${clean_readDir}/${ID}_1.fastq.gz ${clean_readDir}/${ID}_2.fastq.gz)

    # Complete_tag+=(${brackenDir}/sr/${ID}.bracken.out)

    # Complete_tag+=(${shortbredDir}/${ID}.shortbred.tsv)

    # If the completion file exists
    all_exist=true

    for file in "${Complete_tag[@]}"; do
      # echo $file
      if [[ ! -e "$file" ]]; then
        all_exist=false
        break
      fi
    done

    if $all_exist; then
        echo ${ID}
    fi
done


PRE_QC_JOB="Submitted batch job 10501532"
DEDUP_JOB=COMPLETE
TRIM_JOB="Submitted batch job 10501532"
DECONTAM_JOB=COMPLETE
KRAKEN_JOB=COMPLETE
ShortBRED_JOB=COMPLETE

DEPENDENT_JOBs=(${PRE_QC_JOB##* } ${DEDUP_JOB##* } ${TRIM_JOB##* } ${DECONTAM_JOB##* } ${KRAKEN_JOB##* } ${ShortBRED_JOB##* })

ready_to_run=true
dependencies=""
for job in "${DEPENDENT_JOBs[@]}"; do
    ready_to_run=false
    if [[ $job != "COMPLETE" ]]; then
        if [[ -z "$dependencies" ]]; then
            dependencies="${job}"
        else
            dependencies+=",${job}"
        fi
    fi
done

if $ready_to_run; then
    #statements
fi