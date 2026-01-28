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
    H1 "Description: 5.2_NT_amr_bins.sh"
        echo -e "This script does the following:"
        echo -e "1. Combines all MAGs into a single file"
        echo -e "2. Performs AMR gene annotation of .fna files using AMRFinderPlus"
        echo -e "3. Performs AMR gene annotation of .fna files using RGI"

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
        inputType=bin

        H2 "Input"
            if [[ $inputType == "asm" ]]; then
                inputFile=${evaluationDir}/${ID}_final_assembly.fasta
                echo -e $inputFile
            fi

            if [[ $inputType == "bin" ]]; then
                echo -e "${reassemDir}/${ID}/reassembled_bins"
            fi


        H2 "Output"
            if [[ ! -d ${moduleDir}/AMR ]]; then mkdir -p ${moduleDir}/AMR ; fi
            if [[ ! -d ${moduleDir}/RGI ]]; then mkdir -p ${moduleDir}/RGI ; fi

            amrfinder_output=${moduleDir}/AMR/${ID}.amrfinder.txt
            rgi_output=${moduleDir}/RGI/${ID}.rgi.txt

            echo -e "${amrfinder_output}"
            echo -e "${rgi_output}"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source $miniforge_init
    # source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    # conda activate metawrap-env

# Run functions ------------------------------------------------------------------------------------------------------------
all_output_exists=$(test_for_output "${amrfinder_output}" "${rgi_output}")
if ! all_output_exists; then
    
    # Combine all reassembled bins
        if [[ $inputType == "bin" ]]; then
            inputFile=${moduleDir}/${ID}_all_reassembled_bins.fa

            # If ${moduleDir}/${ID}_all_reassembled_bins.fa doesn't exist
            all_bins_exists=$(test_for_output "${inputFile}")
            if  ! $all_bins_exists; then
                H1 "Combine all reassembled bins"
                cat "${reassemDir}/${ID}/reassembled_bins"/*.fa > "${moduleDir}/${ID}_all_reassembled_bins.fa"
            fi
            substep_completion "${inputFile}"
            
        fi


    # AMRFinderPlus
        func="AMRFinderPlus"
        H1 "$func"

        amr_output_exists=$(test_for_output "${amrfinder_output}")

        if $amr_output_exists; then
            comment "Output file already found. Skipping this command..."
        else
            start=$SECONDS
            #------------

            module load blast/2.11.0+
            module load hmmer/3.3.2
            conda activate $AMRFINDER_ENV

            amrfinder \
                --threads $SLURM_CPUS_PER_TASK \
                -n $inputFile \
                -o ${amrfinder_output}

            conda deactivate
            module unload blast/2.11.0+
            module unload hmmer/3.3.2

            #------------
            end=$SECONDS; duration=$(( end-start ))
           
        fi
         step_completion "${amrfinder_output}"


    # RGI
        func="RGI"
        H1 "$func"
        rgi_output_exists=$(test_for_output "${rgi_output}")

        if $rgi_output_exists; then
            comment "Output file already found. Skipping this command..."
        else
            start=$SECONDS
            #------------

            module load diamond/2.0.9
            conda activate $RGI_ENV

            rgi main \
                -i $inputFile \
                -t contig \
                -o ${moduleDir}/RGI/${ID}.rgi \
                --clean \
                -n $SLURM_CPUS_PER_TASK

            conda deactivate
            module unload diamond/2.0.9

            #------------
            end=$SECONDS; duration=$(( end-start ))

        fi
        step_completion "${rgi_output}"


fi


# Completion
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then 
        touch ${moduleDir}/COMPLETE/$ID

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm $int_file
        done
    else
        error "Not all outputs were created."
    fi



H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
