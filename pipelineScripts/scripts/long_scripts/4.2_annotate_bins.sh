#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Performs functional annotation of bin sets with Prokka using metaWRAP's annotate_bins module:
    #     1. Shortens contig names to run Prokka (shorten_contig_names.py).
    #     2. Runs Prokka.
    #     3. Pulls out functional annotations for each bin.
    #     4. Filters the translated genes.
    #     5. Filters the untranslated genes.

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

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
    source $print_functions


# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 4.2_annotate_bins.sh"
        echo -e "This script does the following:"
        echo -e "1. Performs functional annotation of bin sets with Prokka using metaWRAP's annotate_bins module:"
        echo -e "    1. Shortens contig names to run Prokka (shorten_contig_names.py)."
        echo -e "    2. Runs Prokka."
        echo -e "    3. Pulls out functional annotations for each bin."
        echo -e "    4. Filters the translated genes."
        echo -e "    5. Filters the untranslated genes."

    H1 "Job Context"
        OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))
        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

        echo
        print "----- Resources Requested -----"
        print "Nodes:            $SLURM_NNODES"
        print "Cores / node:     $SLURM_CPUS_PER_TASK"
        print "Total memory:     $Total_Gb Gb"
        print "Wall-clock time:  $JobTime"
        print "-------------------------------"

    H1 "Variables"
        comment "SampleID (ID): ${ID}"

        H2 "Input"
            BINS=${reassemDir}/${ID}/reassembled_bins; echo -e "$BINS"

        H2 "Output"
            out=${moduleDir}/$ID
            if [[ ! -d ${out} ]]; then mkdir -p $out; fi
            echo -e "${out}/bin_funct_annotations/"

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

# Run script ---------------------------------------------------------------------------------------------------------------
H1 "Annotate bins"
    start=$SECONDS

    metawrap annotate_bins -o $out -t $SLURM_CPUS_PER_TASK -b $BINS

    conda deactivate
    module unload anaconda3/2023.09

    if [[ $(ls ${out}/bin_funct_annotations | wc -l) -eq 0 ]]; then error "Something went wrong with the annotation."; fi

# Completion status
    if [[ -d ${out}/bin_funct_annotations ]]; then
        touch ${moduleDir}/COMPLETE/${ID}

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm -rf $int_file
        done
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
