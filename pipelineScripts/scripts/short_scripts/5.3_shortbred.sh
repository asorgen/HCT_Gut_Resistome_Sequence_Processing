#!/bin/bash

# Description
    #-----------------------------------------------
    # This script is used to perform AMR gene annotations using ShortBRED.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on June 26, 2025

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

# Config files -------------------------------------------------------------------------------------------------------------
    source $pipelineConfig
    source $config_file
    source $bashrc
    source $bash_profile

# Set function for output comments -----------------------------------------------------------------------------------------
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }
    test_for_output() {
        File_list=("$@")  # capture all arguments as an array
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done
        echo "$all_complete"
    }
    step_completion() {
        File_list=("$@")  # capture all arguments as an array
        
        all_complete=true
        for file in "${File_list[@]}"; do 
            if [[ ! -e "$file" ]]; then 
                all_complete=false
                break
            fi
        done

        if $all_complete; then 
            comment "SUCCESS: $STEP Complete"
        else
            error "[ $STEP ERROR! ] - Exiting..."
        fi
        Complete_tag+=("${File_list[@]}")
    }
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
        comment "This script is used to perform AMR gene annotations using ShortBRED."

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
        echo -e "SampleID (ID): ${ID}"

        H2 "Input"
        FASTQ=${clean_readDir}/${ID}
        echo -e "${FASTQ}_1.fastq.gz"
        echo -e "${FASTQ}_2.fastq.gz"
        OUTPUT_MARKERS=$SB_MARKERS
        echo -e $OUTPUT_MARKERS
        TEMP_DIR="${moduleDir}/${ID}_tmp"
        out=$moduleDir

        H2 "Output"
        OUTPUT_FILE="${out}/${ID}.shortbred.tsv"
        echo -e "$OUTPUT_FILE"

        if [[ ! -d ${out}/COMPLETE ]]; then mkdir -p $out/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load shortbred/0.9.5

# Run functions ------------------------------------------------------------------------------------------------------------
    H2 "ShortBRED Quantify"

    shortbred_quantify.py \
        --markers $OUTPUT_MARKERS \
        --wgs ${FASTQ}_1.fastq.gz ${FASTQ}_2.fastq.gz \
        --results $OUTPUT_FILE \
        --tmp $TEMP_DIR
    Complete_tag+=($OUTPUT_FILE)

    if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi

    rm -r "$TEMP_DIR"

    H1 "Completion"
        
        output_exists=$(test_for_output "${Complete_tag[@]}")
        if $output_exists; then touch ${out}/COMPLETE/$ID; fi

# Unload environments ------------------------------------------------------------------------------------------------------
    module unload shortbred/0.9.5

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"

