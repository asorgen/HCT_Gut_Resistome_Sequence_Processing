#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Performs taxonomic relative abundance profiling with MetaPhlAn4.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 7, 2024

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
    substep_completion() {
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
        Intermediate_files+=("${File_list[@]}")
    }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
        echo -e "This script does the following:"
        echo -e "1. Performs taxonomic relative abundance profiling with MetaPhlAn4."

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
            R1=${clean_readDir}/${ID}_1.fastq.gz; echo -e "${R1}"
            R2=${clean_readDir}/${ID}_2.fastq.gz; echo -e "${R2}"
            if [[ -s ${evaluationDir}/${ID}_final_assembly.fasta ]]; then ASM=${evaluationDir}/${ID}_final_assembly.fasta; echo -e "$ASM"; fi
            
        H2 "Output"
            out=${moduleDir}
            echo -e "${out}/${ID}_rel_abun_w_read_stats.txt"
            # echo -e "${out}/${ID}_clade_profile.txt"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p $moduleDir/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load metaphlan/4.2.2
    
# Run functions ------------------------------------------------------------------------------------------------------------

STEP="MetaPhlAn4 Relative Abundance Profiling"
    
    H1 "$STEP"

      RA_files=("${out}/${ID}_rel_abun_w_read_stats.txt")
      output_exists=$(test_for_output "${RA_files[@]}")

      if ! $output_exists; then

        CMD="metaphlan $R1,$R2 --input_type fastq -t rel_ab_w_read_stats --force --db_dir $METAPHLAN4_DB --nproc $SLURM_CPUS_PER_TASK --output_file ${out}/${ID}_rel_abun_w_read_stats.txt --mapout ${out}/${ID}_mapping.txt --ignore_eukaryotes --ignore_archaea"
        comment "$CMD"
        $CMD
        if [[ $? -ne 0 ]]; then error "Something went wrong! Exiting..."; fi

      fi
      step_completion "${RA_files[@]}"


# STEP="MetaPhlAn4 Clade Profiling"
    
#     H1 "$STEP"

#       CP_files=("${out}/${ID}_clade_profile.txt" "${out}/${ID}_mapping.txt")
#       output_exists=$(test_for_output "${CP_files[@]}")

#       if ! $output_exists; then

#         CMD="metaphlan $R1,$R2 --input_type fastq -t clade_profiles --force --db_dir $METAPHLAN4_DB --nproc $SLURM_CPUS_PER_TASK --output_file ${out}/${ID}_clade_profile.txt --mapout ${out}/${ID}_mapping.txt --ignore_eukaryotes --ignore_archaea"
#         echo -e $CMD
#         $CMD
#         if [[ $? -ne 0 ]]; then error "Something went wrong! Exiting..."; fi

#       fi
#       step_completion "${CP_files[@]}"


H1 "Completion"
    
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then 
        touch ${moduleDir}/COMPLETE/$ID

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm $int_file
        done
    fi

# Unload environments ------------------------------------------------------------------------------------------------------
    module unload metaphlan/4.2.2

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"



