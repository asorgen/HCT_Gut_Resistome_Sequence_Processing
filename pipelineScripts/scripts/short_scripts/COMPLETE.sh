#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script is used to all intermediate files for finished pipeline samples.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Oct 13, 2025

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

# Slurm Resource Options ---------------------------------------------------------------------------------------------------

    # Job partition (--partition=<partition_names>; -p <partition_names>; SBATCH_PARTITION) | Options: Orion, Nebula, Pisces
    # Job name (--job-name=<name>; -J <name>; SBATCH_JOB_NAME)
    # Path to file storing text output. (--output=<filename_pattern>; -o <name>; SBATCH_OUTPUT)
    # Node count required for the job (--nodes=<count>; -N <count>)
    # Request that ntasks be invoked on each node. (--ntasks-per-node=<ntasks>)
    # Memory required per node (--mem=<MB>[units]; SLURM_MEM_PER_NODE)
    # Notify user by email when certain event types occur. (--mail-type=<type>) | Options: NONE, BEGIN, END, FAIL, REQUEUE, ALL
    # User to receive email notification of state changes as defined by --mail-type. (--mail-user=<user>)
    # Maximum allowed runtime of job (--time=<time>; -t <time>; SBATCH_TIMELIMIT)


source ${HOME}/.bashrc
source $config_file
source $pipelineConfig

# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }
    test_for_output() {
        File_list=("$@")  # capture all arguments as an array
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done
        echo "$all_complete"
    }

    step_completion() {
        File_list=("$@")  # capture all arguments as an array
        
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done

        if $all_complete; then 
            comment "SUCCESS: $STEP Complete"
        else
            error "[ $STEP ERROR! ] - Exiting..."
        fi
        Complete_tag+=("${File_list[@]}")
    }

H1 "Usage"
    comment "This script is used to all intermediate files for finished pipeline samples."

H1 "Job Context"
    export SLURM_NTASKS
    comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
    comment "Running on host: `hostname`"

    Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))

    JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

    echo 
    comment "----- Resources Requested -----"
    comment "Nodes:            $SLURM_NNODES"
    comment "Cores / node:     $SLURM_NTASKS"
    comment "Total memory:     $Total_Gb Gb"
    comment "Wall-clock time:  $JobTime"
    comment "-------------------------------"

H1 "Variables"
    echo -e "SampleID (ID): ${ID}"

    H2 "Output"
        out=${datasetDir}/COMPLETE
        mkdir -p ${out}
        echo ${out}/${ID}


H3 "[ Start ]"
/bin/date
SECONDS=0
start=$SECONDS
Complete_tag=()
Intermediate_files=()

H1 "Double check outputs exist"

# 0.0_raw_reads
Intermediate_files+=(${raw_readDir}/${ID}_1.fastq)
Intermediate_files+=(${raw_readDir}/${ID}_2.fastq)

# 0.1_pre_qc

# 0.2_deduplication
Intermediate_files+=(${dedup_Dir}/${ID}_deduped_R1.fastq.gz)
Intermediate_files+=(${dedup_Dir}/${ID}_deduped_R2.fastq.gz)

# 0.3_sequence_trim
Intermediate_files+=(${trimmed_Dir}/${ID}_trimmed_1.fastq.gz)
Intermediate_files+=(${trimmed_Dir}/${ID}_trimmed_2.fastq.gz)

# 0.4_host_decontamination
STEP="Sequence Processing"
Final_Fastqs=(${clean_readDir}/${ID}_1.fastq.gz ${clean_readDir}/${ID}_2.fastq.gz)
step_completion "${Final_Fastqs[@]}"

# 2.1_kraken2
STEP="Kraken2"
Intermediate_files+=(${krakenDir}/${ID}/${ID}.krak2 ${krakenDir}/${ID}/${ID}.kraken2)

# 2.2_bracken
STEP="Taxonomic Classification"
step_completion "${brackenDir}/sr/${ID}.bracken.out"

# 5.3_shortbred
if $run_shortbred; then
    STEP="ShortBRED"
    step_completion "${shortbredDir}/${ID}.shortbred.tsv"
fi

# 5.4_rgi_bwt
if $run_rgi_bwt; then
    STEP="RGI BWT"
    step_completion "${rgi_bwt_dir}/kma_output/${ID}.rgi_kma.txt"
fi

H1 "Completion"
    
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then 
        CMD="touch ${out}/${ID}"; echo $CMD
        $CMD

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            if [[ -e $int_file ]]; then
                CMD="rm $int_file"; echo $CMD
                $CMD
            fi
            
        done
    else
        error "Pipeline failed"
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
