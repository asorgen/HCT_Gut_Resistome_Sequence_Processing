#!/bin/bash

# Description
    #-----------------------------------------------
    # This script is used to assemble hybrid-read metagenome samples.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 7, 2024

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1
    #-----------------------------------------------

# Slurm Resource Options
    #-----------------------------------------------

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

#SBATCH --mail-user=asorgen@uncc.edu

source /users/asorgen/.bashrc
# config_file=$(which config-metawrap)
# source $config_file

# Set function for output comments
H1 () { print_header.py "$1" "H1"; }
H2 () { print_header.py "$1" "H2"; }
H3 () { print_header.py "$1" "H3"; }
comment () { print_header.py "$1" "#"; }
error () { echo $1; exit 1; }

H1 "Usage"
    comment "This script is used to assemble Oxford Nanopore long read metagenome samples."

H1 "Variables"
    comment "SampleID (ID): ${ID}"
    H2 "Input"
    echo -e "../SHORT/${clean_readDir}/${ID}_1.fastq"
    echo -e "../SHORT/${clean_readDir}/${ID}_2.fastq"
    echo -e "../LONG/${clean_readDir}/${ID}_ont.fastq"
    H2 "Output"
    echo -e "Assembly files will be deposited to ${assemblyDir}/${ID}"
    mkdir -p ${assemblyDir}/${ID}
    

H1 "Job Context"
    OMP_NUM_THREADS=$SLURM_NTASKS
    comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
    comment "Running on host: `hostname`"

    Total_Gb=$(( SLURM_MEM_PER_NODE / 1000 ))

    JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

    echo 
    comment "----- Resources Requested -----"
    comment "Nodes:            $SLURM_NNODES"
    comment "Cores / node:     $SLURM_NTASKS"
    comment "Total memory:     $Total_Gb Gb"
    comment "Wall-clock time:  $JobTime"
    comment "-------------------------------"

H2 "[ Start ]"
/bin/date
SECONDS=0

export PATH="/users/asorgen/PROGRAMS/metaWRAP/bin/metawrap-scripts/:$PATH"

func="OPERA-MS Assembly"
    H1 "$func"

    export PERL5LIB="$HOME/perl5/lib/perl5"
    module load R/4.2.2

    perl /users/asorgen/PROGRAMS/OPERA-MS/OPERA-MS.pl \
        --num-processors 32 \
        --short-read1 ${clean_readDir}/${ID}_1.fastq \
        --short-read2 ${clean_readDir}/${ID}_2.fastq \
        --long-read ${clean_readDir}/${ID}_ont.fastq \
        --out-dir ${assemblyDir}/${ID}
    if [[ $? -ne 0 ]]; then error "Something went wrong with the assembly. Exiting"; fi

    # Default metaFlye assembly file: contigs.polished.fasta
    mv ${assemblyDir}/${ID}/contigs.polished.fasta ${assemblyDir}/${ID}/assembly.fasta



H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"


