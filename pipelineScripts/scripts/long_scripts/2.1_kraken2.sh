#!/bin/bash

#-----------------------------------------------
# This script is used to perform taxonomic classifications on reads and assemblies of Illumina metagenome sequences.

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

#SBATCH --mail-user=${email}

source ${HOME}/.bashrc
metawrap_config=$(which config-metawrap)
source $metawrap_config

# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }

H1 "Usage"
    comment "This script is used to perform taxonomic classifications on reads and assemblies of Illumina metagenome sequences."

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

H1 "Variables"
    echo -e "SampleID (ID): ${ID}"
    H2 "Input"
        ONT=${clean_readDir}/${ID}_ont.fastq
        R1=${clean_readDir}/${ID}_1.fastq
        R2=${clean_readDir}/${ID}_2.fastq
        ASM=${evaluationDir}/${ID}_final_assembly.fasta

        if [[ "$readType" == "long_read" ]]; then
            echo -e "${ONT}"
            echo -e "$ASM"
        fi

        if [[ "$readType" == "hybrid" ]]; then
            echo -e "$R1"
            echo -e "$R2"
            echo -e "$ONT"
            echo -e "$ASM"
        fi

    H2 "Output"
        k_out=${moduleDir}/${ID}
        b_out=bracken
        echo -e "Kraken2 and Krona output will be deposited to ${k_out}/"
        echo -e "Bracken will be deposited to ${b_out}/"



H2 "[ Start ]"
/bin/date
SECONDS=0

H1 "Kraken2"
    start=$SECONDS

    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda init
    conda activate metawrap-env

    if [[ "$readType" == "hybrid" ]]; then
        sr_out=${k_out}/${ID}.krak2
        sr_rep=${k_out}/${ID}.kreport

        
        if [ ! -s $sr_out ] || [ ! -s $sr_rep ]; then
            H2 "Processing short reads"
            output=$sr_out; report=$sr_rep
            CMD="kraken2 --use-names --db ${KRAKEN2_DB} --paired --threads $SLURM_NTASKS --output $output --report $report $R1 $R2"
            echo $CMD
            $CMD
        fi

        if [[ $? -ne 0 ]] || [[ ! -s $sr_out ]] || [[ ! -s $sr_rep ]]; then error "Something went wrong with Kraken2 while processing the short-read fastqs. Exiting..."; fi
    fi

    ont_out=${k_out}/${ID}_ONT.krak2
    ont_rep=${k_out}/${ID}_ONT.kreport

    asm_out=${k_out}/${ID}_assembly.krak2
    asm_rep=${k_out}/${ID}_assembly.kreport

    
    if [ ! -s $ont_out ] || [ ! -s $ont_rep ]; then
        H2 "Processing ONT reads"
        output=$ont_out; report=$ont_rep; input=$ONT
        CMD="kraken2 --use-names --db ${KRAKEN2_DB} --threads $SLURM_NTASKS --output $output --report $report $input"
        echo $CMD
        $CMD
    fi

    if [[ $? -ne 0 ]] || [[ ! -s $ont_out ]] || [[ ! -s $ont_rep ]]; then error "Something went wrong with Kraken2 while processing the ONT fastq. Exiting..."; fi

    if [ ! -s $asm_out ] || [ ! -s $asm_rep ]; then
        H2 "Processing assembly"
        output=$asm_out; report=$asm_rep; input=$ASM
        CMD="kraken2 --use-names --db ${KRAKEN2_DB} --threads $SLURM_NTASKS --output $output --report $report $input"
        echo $CMD
        $CMD
    fi

    if [[ $? -ne 0 ]] || [[ ! -s $asm_out ]] || [[ ! -s $asm_rep ]]; then error "Something went wrong with Kraken2 while processing the assembly fasta. Exiting..."; fi

    
H1 "Kraken2 Translate"
    for file in ${k_out}/*.krak2; do
        comment "Translating $file"
        ${SOFT}/kraken2_translate.py ${KRAKEN2_DB} $file ${file%.*}.kraken2
        if [[ $? -ne 0 ]]; then error "Something went wrong with kraken2-translate. Exiting..."; fi
    done

H1 "Kronogram"
    for file in ${k_out}/*.kraken2; do
        ${SOFT}/kraken_to_krona.py $file > ${file%.*}.krona
        if [[ $? -ne 0 ]]; then error "Something went wrong while making the krona file. Exiting..."; fi
    done

    ktImportText -o ${k_out}/kronagram.html ${k_out}/*krona
    if [[ ! -s ${k_out}/kronagram.html ]]; then error "Something went wrong while running KronaTools to make kronagram. Exiting..."; fi

    conda deactivate
    module unload anaconda3/2023.09


H1 "Bracken"

    module load anaconda3/2023.09

    for file in ${k_out}/*.kreport; do
        python3 ${HOME}/PROGRAMS/Bracken-2.7/src/est_abundance.py -i $file -k ${KRAKEN2_DB}/database150mers.kmer_distrib --level S -o ${file%.*}.bracken.out
        if [[ $? -ne 0 ]]; then error "Something went wrong while running Bracken. Exiting..."; fi
    done

    mkdir -p ${b_out}/assembly
    mv ${k_out}/*_assembly.bracken.out ${b_out}/assembly
    
    mkdir -p ${b_out}/ont
    mv ${k_out}/*_ONT.bracken.out ${b_out}/ont

    if [[ "$readType" == "hybrid" ]]; then
        mkdir -p ${b_out}/sr
        mv ${k_out}/*.bracken.out ${b_out}/sr
    fi


    module unload anaconda3/2023.09

# Completion status
    if [[ -s bracken/assembly/${ID}_assembly.bracken.out && -s bracken/ont/${ID}_ONT.bracken.out ]]; then
        
        if [[ "$readType" == "hybrid" && -s bracken/sr/${ID}.bracken.out ]]; then
            touch ${moduleDir}/${ID}/COMPLETE
        fi

        if [[ "$readType" == "long_read" ]]; then
            touch ${moduleDir}/${ID}/COMPLETE
        fi
        
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"



