#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script removes human host reads from trimmed ONT reads using minimap2 and samtools:
    # 1. Aligns trimmed reads to GRCh38 with minimap2 in map-ont mode
    # 2. Extracts unmapped (non-human) reads with samtools
    # 3. Outputs clean reads as ${clean_readDir}/${ID}_ont.fastq

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

    # ONT-specific notes:
    #   - Replaces Bowtie2 + GRCh38 from the short-read pipeline
    #   - minimap2 -ax map-ont is optimized for Oxford Nanopore reads
    #   - samtools -f 4 retains only unmapped (non-human) reads
    #   - Requires GRCh38_FASTA to be set in the config (path to GRCh38 reference FASTA)
    #   - Input: ${trimmed_Dir}/${ID}_trimmed.fastq.gz
    #   - Output: ${clean_readDir}/${ID}_ont.fastq (uncompressed, for downstream compatibility)

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
    H1 "Description: 0.4_host_decontamination.sh (ONT)"
        echo -e "This script removes human host reads from trimmed ONT reads:"
        echo -e "1. Aligns reads to GRCh38 with minimap2 (map-ont mode)"
        echo -e "2. Extracts unmapped (non-human) reads with samtools"
        echo -e "3. Outputs: ${clean_readDir}/${ID}_ont.fastq"

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
            echo -e "${trimmed_Dir}/${ID}_trimmed.fastq.gz"
            echo -e "GRCh38 reference: ${GRCh38_FASTA}"
        H2 "Output"
            if [[ ! -d ${clean_readDir} ]]; then mkdir -p ${clean_readDir}; fi
            echo -e "${clean_readDir}/${ID}_ont.fastq"
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi
            mkdir -p ${moduleDir}/QC_report

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env

    module load samtools

# Run functions ------------------------------------------------------------------------------------------------------------

STEP="Host read removal (minimap2 map-ont + samtools)"
    H1 "$STEP"

    outputFile=${clean_readDir}/${ID}_ont.fastq
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        H3 "Aligning to GRCh38 with minimap2 (map-ont) and extracting unmapped reads"
        minimap2 \
            -ax map-ont \
            -t $SLURM_CPUS_PER_TASK \
            ${GRCh38_FASTA} \
            ${trimmed_Dir}/${ID}_trimmed.fastq.gz \
        | samtools view -bS -f 4 - \
        | samtools sort -@ $SLURM_CPUS_PER_TASK - \
        | samtools fastq - > ${clean_readDir}/${ID}_ont.fastq

        if [[ $? -ne 0 ]]; then error "Host decontamination failed for ${ID}. Exiting..."; fi

        H3 "Read count summary"
        input_count=$(( $(zcat ${trimmed_Dir}/${ID}_trimmed.fastq.gz | wc -l) / 4 ))
        clean_count=$(( $(cat ${clean_readDir}/${ID}_ont.fastq | wc -l) / 4 ))
        host_count=$(( input_count - clean_count ))
        comment "Trimmed reads input:    ${input_count}"
        comment "Host reads removed:     ${host_count}"
        comment "Clean reads remaining:  ${clean_count}"

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$outputFile" ]]; then
            H2 "Host Decontamination Complete"
            comment "$(elapsed_time "$duration")"
        fi
    fi
    step_completion "${outputFile}"


STEP="Post-decontamination QC (FastQC)"
    H1 "$STEP"

    qc_out=${moduleDir}/QC_report/${ID}_fastqc.zip
    if [[ -s "$qc_out" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        fastqc -q -t $SLURM_CPUS_PER_TASK \
            -o ${moduleDir}/QC_report \
            -f fastq ${clean_readDir}/${ID}_ont.fastq

        if [[ $? -ne 0 ]]; then echo "Warning: FastQC post-decontam failed for ${ID}. Continuing..."; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$qc_out" ]]; then
            H2 "Post-QC Complete"
            comment "$(elapsed_time "$duration")"
        fi
    fi

    conda deactivate
    module unload samtools
    module unload anaconda3/2023.09

# Completion
    module_completion

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
