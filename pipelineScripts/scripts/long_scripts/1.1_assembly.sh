#!/bin/bash

#-----------------------------------------------
# This script is used to assemble Oxford Nanopore long read metagenome samples.

# It is important to note that it has been designed for a specific working directory. 
# Therefore, the reproduction of the results will require small modifications of the script 
# or the adaptation of your working directory.

# Created on Nov 7, 2024

# @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

# Version: 1
#-----------------------------------------------

#-----------------------------------------------
# Slurm Resource Options

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

##SBATCH --mail-user=${email}

source $pipelineConfig
source $config_file
source $bashrc
source $bash_profile

# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }

H1 "Usage"
    comment "This script is used to assemble Oxford Nanopore long read metagenome samples."

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
    echo -e "${clean_readDir}/${ID}_ont.fastq"
    H2 "Output"
    echo -e "Assembly files will be deposited to ${moduleDir}/${ID}"
    mkdir -p ${moduleDir}/${ID}
    

H2 "[ Start ]"
/bin/date
SECONDS=0



func="metaFlye Assembly"
    H1 "$func"
    
    
    # Force Lustre metadata refresh before checking for existing assembly
    ls "${moduleDir}/${ID}/" > /dev/null 2>&1 || true

    if [[ ! -s "${moduleDir}/${ID}/assembly.fasta" ]]; then
        module load flye
        if [[ -f "${moduleDir}/${ID}/params.json" ]]; then
            flye --nano-hq ${clean_readDir}/${ID}_ont.fastq --out-dir ${moduleDir}/${ID} --meta -t $SLURM_CPUS_PER_TASK --resume
        else
            flye --nano-hq ${clean_readDir}/${ID}_ont.fastq --out-dir ${moduleDir}/${ID} --meta -t $SLURM_CPUS_PER_TASK
        fi
        
        if [[ $? -ne 0 ]]; then echo "Something went wrong with ${func}. Exiting"; exit 1; fi

        # Default metaFlye assembly file: assembly.fasta

        module unload flye
    fi

# func="Renaming contigs"
#     H1 "$func"

#     module load anaconda3/2023.09 
#     cp ${moduleDir}/${ID}/assembly.fasta ${moduleDir}/${ID}/orig_assembly.fasta
#     scripts/helper_scripts/rename_flye_assembly.py ${moduleDir}/${ID}/orig_assembly.fasta ${moduleDir}/${ID}/assembly_info.txt > ${moduleDir}/${ID}/renamed_contigs.fasta

# Completion status
    if [[ -s ${moduleDir}/${ID}/assembly.fasta ]]; then
        touch ${moduleDir}/COMPLETE/${ID}
    fi


H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"


