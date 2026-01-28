#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script is used to reassemble bins according to a "permissive" and a "strict" algorithm. 
    # Only bins that have been improved through reassembly will be altered in the final set.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 10, 2025

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
    pFunc () { echo $1; echo; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
        comment "This script is used to reassemble bins according to a \"permissive\" and a \"strict\" algorithm. Only bins that have been improved through reassembly will be altered in the final set."

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
            reads_1=${clean_readDir}/${ID}_1.fastq.gz; echo -e "$reads_1"
            reads_2=${clean_readDir}/${ID}_2.fastq.gz; echo -e "$reads_2"

            if [[ "$readType" == "hybrid" ]]; then
                reads_ONT=${clean_readDir}/${ID}_ont.fastq.gz; echo -e "$reads_ONT"
            fi

            BINS_refined=${refinedbinDir}/${ID}/metawrap_${min_completion}_${max_contam}_bins; echo -e "$BINS_refined"

        H2 "Output"
            out=${moduleDir}/$ID
            if [[ ! -d ${out} ]]; then mkdir -p $out; fi
            BINS_reassembled=${out}/reassembled_bins; echo -e "$BINS_reassembled"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env

STEP="Unzip and rename the sequence files"
    R1=${out}/${ID}_1.fastq
    R2=${out}/${ID}_2.fastq
    ONT=${out}/${ID}_ont.fastq

    Unzip_files=(${R1} ${R2})
    if [[ "$readType" == "hybrid" ]]; then Unzip_files+=(${ONT}); fi
    
    output_exists=$(test_for_output "${Unzip_files[@]}")
    if ! "$output_exists"; then
        H1 "$STEP"
        echo -e "gunzip -c ${reads_1} > $R1"
        gunzip -c ${reads_1} > $R1

        echo -e "gunzip -c ${reads_2} > $R2"
        gunzip -c ${reads_2} > $R2
        if [[ "$readType" == "hybrid" ]]; then gunzip -c ${reads_ONT} > $ONT; fi
    fi 
    substep_completion "${Unzip_files[@]}"

H1 "Reassemble bins"
    start=$SECONDS

    if [[ "$readType" == "hybrid" ]]; then
        metawrap reassemble_bins \
            -o $out \
            -1 $R1 \
            -2 $R2 \
            --nanopore $ONT \
            -t $SLURM_CPUS_PER_TASK \
            -m $Total_Gb \
            -c $min_completion \
            -x $max_contam \
            -b $BINS_refined
    else
        metawrap reassemble_bins \
            -o $out \
            -1 $R1 \
            -2 $R2 \
            -t $SLURM_CPUS_PER_TASK \
            -m $Total_Gb \
            -c $min_completion \
            -x $max_contam \
            -b $BINS_refined
    fi

    if [[ ! -s ${moduleDir}/${ID}/reassembled_bins.stats ]]; then error "Something went wrong with bin reassembly. Exiting..."; fi

    Intermediate_files+=(${out}/original_bins)
    Intermediate_files+=(${out}/work_files)
    Intermediate_files+=(${out}/reassembled_bins.checkm)

    conda deactivate
    module unload anaconda3/2023.09

    # File containing the data
    file="${out}/reassembled_bins.stats"

    # Count bins ending in .strict
    strict_count=$(grep -c '\.strict' "$file")

    # Count bins ending in .permissive
    permissive_count=$(grep -c '\.permissive' "$file")

    # Count bins ending in .orig
    orig_count=$(grep -c '\.orig' "$file")

    H2 "Reassembly Results"
    # Output the results
    comment "Bins improved through strict reassembly: $strict_count"
    comment "Bins improved through permissive reassembly: $permissive_count"
    comment "Bins that could not be improved: $orig_count"
    
# Completion status
    if [[ -s ${moduleDir}/${ID}/reassembled_bins.stats ]]; then
        touch ${moduleDir}/COMPLETE/${ID}

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm -rf $int_file
        done

    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"



