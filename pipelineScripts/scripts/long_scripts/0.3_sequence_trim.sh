#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following for Oxford Nanopore Technology (ONT) reads:
    # 1. Adapter trimming with Porechop
    # 2. Quality and length filtering with NanoFilt

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

    # ONT-specific notes:
    #   - Replaces TrimGalore/Cutadapt from the short-read pipeline
    #   - Porechop removes ONT-specific adapters (barcodes, ligation adapters)
    #   - NanoFilt filters by minimum quality score and minimum read length
    #   - Input: ${raw_readDir}/${ID}.fastq.gz (single-end, gzipped)
    #   - Output: ${trimmed_Dir}/${ID}_trimmed.fastq.gz

# Slurm Resource Options ---------------------------------------------------------------------------------------------------

    # Job partition (--partition=<partition_names>; -p <partition_names>; SBATCH_PARTITION) | Options: Orion, Nebula, Pisces
    # Job name (--job-name=<name>; -J <name>; SBATCH_JOB_NAME)
    # Path to file storing text output. (--output=<filename_pattern>; -o <name>; SBATCH_OUTPUT)
    # Node count required for the job (--nodes=<count>; -N <count>)
    # Request that cpus per task (--cpus-per-task=<ncpus>)
    # Memory required per node (--mem=<MB>[units]; SLURM_MEM_PER_NODE)
    # Notify user by email when certain event types occur. (--mail-type=<type>) | Options: NONE, BEGIN, END, FAIL, REQUEUE, ALL
    # User to receive email notification of state changes as defined by --mail-type. (--mail-user=<user>)
    # Maximum allowed runtime of job (--time=<time>; -t <time>; SBATCH_TIMELIMIT)

# Config files -------------------------------------------------------------------------------------------------------------
    source $pipelineConfig
    source $config_file
    source $bashrc
    source $bash_profile
    source $module_functions

# Set function for output comments -----------------------------------------------------------------------------------------
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 0.3_sequence_trim.sh (ONT)"
        echo -e "This script does the following:"
        echo -e "1. Adapter trimming with Porechop"
        echo -e "2. Quality and length filtering with NanoFilt (Q>=${nanofilt_min_quality}, len>=${nanofilt_min_length}bp)"

    H1 "Job Context"
        OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))
        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

        echo
        comment "----- Resources Requested -----"
        comment "Nodes:            $SLURM_NNODES"
        comment "Cores / node:     $SLURM_CPUS_PER_TASK"
        comment "Total memory:     $Total_Gb Gb"
        comment "Wall-clock time:  $JobTime"
        comment "-------------------------------"

    H1 "Variables"
        comment "SampleID (ID): ${ID}"
        H2 "Input"
            echo -e "${raw_readDir}/${ID}.fastq.gz"
        H2 "Output"
            if [[ ! -d ${trimmed_Dir} ]]; then mkdir -p ${trimmed_Dir}; fi
            echo -e "${trimmed_Dir}/${ID}_trimmed.fastq.gz"
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate $ONT_ENV

# Run functions ------------------------------------------------------------------------------------------------------------

STEP="Porechop adapter trimming"
    H1 "$STEP"

    outputFile="${moduleDir}/${ID}_porechop.fastq.gz"
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        porechop_abi \
            -t $SLURM_CPUS_PER_TASK \
            -i ${raw_readDir}/${ID}.fastq.gz \
            -o ${moduleDir}/${ID}_porechop.fastq.gz

        if [[ $? -ne 0 ]]; then error "Porechop failed for ${ID}. Exiting..."; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$outputFile" ]]; then
            H2 "Porechop Complete"
            comment "$(elapsed_time "$duration")"
        fi
    fi
    substep_completion "${outputFile}"


STEP="NanoFilt quality and length filtering"
    H1 "$STEP"

    outputFile="${trimmed_Dir}/${ID}_trimmed.fastq.gz"
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        gunzip -c ${moduleDir}/${ID}_porechop.fastq.gz \
            | NanoFilt -q ${nanofilt_min_quality} -l ${nanofilt_min_length} \
            | gzip > ${trimmed_Dir}/${ID}_trimmed.fastq.gz

        if [[ $? -ne 0 ]]; then error "NanoFilt failed for ${ID}. Exiting..."; fi

        H3 "Read count summary"
        raw_count=$(( $(zcat ${raw_readDir}/${ID}.fastq.gz | wc -l) / 4 ))
        trimmed_count=$(( $(zcat ${trimmed_Dir}/${ID}_trimmed.fastq.gz | wc -l) / 4 ))
        comment "Raw reads:     ${raw_count}"
        comment "Trimmed reads: ${trimmed_count}"

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$outputFile" ]]; then
            H2 "NanoFilt Complete"
            comment "$(elapsed_time "$duration")"
        fi
    fi
    step_completion "${outputFile}"

    # Remove Porechop intermediate
    Intermediate_files+=(${moduleDir}/${ID}_porechop.fastq.gz)

    conda deactivate
    module unload anaconda3/2023.09

# Completion
    module_completion

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
