#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Performs AMR gene annotation of the filtered assembly (.fasta) using:
    #     1. AMRFinderPlus (nucleotide-level search against NCBI AMRFinderPlus database)
    #     2. RGI main (nucleotide-level search against CARD database)

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

    # Input: ${evaluationDir}/${ID}_final_assembly.fasta (nucleotide FASTA)

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
    pFunc () { echo $1; echo; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 5.1_NT_amr_assembly.sh"
        echo -e "This script does the following:"
        echo -e "1. Performs AMR gene annotation of .fna files using AMRFinderPlus"
        echo -e "2. Performs AMR gene annotation of .fna files using RGI"

    H1 "Job Context"
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
            inputFile=${evaluationDir}/${ID}_final_assembly.fasta
            echo -e $inputFile

        H2 "Output"
            if [[ ! -d ${moduleDir}/${ID} ]]; then mkdir -p ${moduleDir}/${ID}; fi
            amr_output=${moduleDir}/${ID}/${ID}.amrfinder.tsv; echo -e "${amr_output}"
            rgi_output=${moduleDir}/$ID/${ID}.rgi.txt; echo -e "${rgi_output}"
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()


STEP="AMRFinderPlus"
    H1 "$STEP"
    source $miniforge_init

    if [[ -s "$amr_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load blast/2.11.0+
        module load hmmer/3.3.2
        conda activate $AMRFINDER_ENV

        amrfinder \
        --threads $SLURM_CPUS_PER_TASK \
        -n ${inputFile} \
        -o ${amr_output}

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        conda deactivate
        module unload blast/2.11.0+
        module unload hmmer/3.3.2

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$amr_output" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${amr_output}"


STEP="RGI"
    H1 "$STEP"
    rgi_output=${moduleDir}/$ID/${ID}.rgi.txt

    if [[ -s "$rgi_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load diamond/2.0.9
        source $RGI_ENV

        rgi load --card_json $CARD_DB/card.json

        rgi main \
        -i ${inputFile} \
        -t contig \
        -o ${moduleDir}/${ID}/${ID}.rgi \
        --clean \
        -n $SLURM_CPUS_PER_TASK

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        conda deactivate
        module unload diamond/2.0.9

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$rgi_output" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${rgi_output}"


# Completion
    module_completion

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
