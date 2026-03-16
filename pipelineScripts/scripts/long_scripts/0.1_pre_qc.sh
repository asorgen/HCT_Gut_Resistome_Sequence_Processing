#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Runs FastQC on raw ONT reads to generate a pre-QC quality report.

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

    # ONT-specific notes:
    #   - Input is a single .fastq.gz file (not paired-end)
    #   - FastQC is applied directly to the compressed raw reads

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
    H1 "Description: 0.1_pre_qc.sh (ONT)"
        echo -e "This script runs FastQC on raw ONT reads to generate a pre-QC quality report."

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
            echo -e "${moduleDir}/${ID}_fastqc.html"
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env

# Run functions ------------------------------------------------------------------------------------------------------------

STEP="Pre-QC Report (FastQC)"
    H1 "$STEP"

    outputFile=${moduleDir}/${ID}_fastqc.zip
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        fastqc -q -t $SLURM_CPUS_PER_TASK \
            --java-mem-size 4096m \
            -o ${moduleDir} \
            -f fastq ${raw_readDir}/${ID}.fastq.gz

        if [[ $? -ne 0 ]]; then error "FastQC failed for ${ID}. Exiting..."; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))
    fi
    step_completion "${outputFile}"

    conda deactivate
    module unload anaconda3/2023.09

# Completion
    module_completion

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
