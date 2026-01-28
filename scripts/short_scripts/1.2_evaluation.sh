#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------

    # This script does the following:
    # 1. Copies the final metaspades assembly
    # 2. Filters out contigs with <1500 bp.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 7, 2024

    # @author: Prisca - UNC Chapel Hill
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
    #-----------------------------------------------

    #SBATCH --mail-user=${email}

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
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }


# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
        echo -e "This script does the following:"
        echo -e "1. Copies the final metaspades assembly."
        echo -e "2. Filters out contigs with <1500 bp."

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
        echo -e "${assemblyDir}/${ID}/assembly.fasta"

        H2 "Output"
        echo -e "Assembly statistics will be deposited to ${moduleDir}/"
        if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi
    

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    start=$SECONDS

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    
# Run functions ------------------------------------------------------------------------------------------------------------

H1 "Copy metaspades assembly and rename it"
    
    if [[ ! -s "${moduleDir}/${ID}_draft_assembly.fasta" ]]; then
        start=$SECONDS
        #------------
        cp ${assemblyDir}/${ID}/assembly.fasta ${moduleDir}/${ID}_draft_assembly.fasta
        if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi
        #------------
        end=$SECONDS; duration=$(( end-start ))
    fi
    step_completion "${moduleDir}/${ID}_draft_assembly.fasta"

H1 "Filter contigs <1500 bp"

    if [[ ! -s "${moduleDir}/${ID}_final_assembly.fasta" ]]; then
        start=$SECONDS
        #------------
        reformat.sh in=${moduleDir}/${ID}_draft_assembly.fasta out=${moduleDir}/${ID}_final_assembly.fasta minlength=1500
        if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi
        #------------
        end=$SECONDS; duration=$(( end-start ))
    fi
    step_completion "${moduleDir}/${ID}_final_assembly.fasta"

    
module_completion
H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
