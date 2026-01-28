#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Removes adapters and trims/filters out low quality reads using TrimGalore.
    # 2. Performs a quality assessment of the deduplicated, trimmed sequence reads using FastQC.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Oct 10, 2025

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
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
        echo -e "This script does the following:"
        echo -e "1. Removes adapters and trims/filters out low quality reads using TrimGalore."
        echo -e "2. Performs a quality assessment of the deduplicated, trimmed sequence reads using FastQC."

    H1 "Job Context"
        export SLURM_NTASKS
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))

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
            R1=${dedup_Dir}/${ID}_deduped_R1.fastq.gz; echo -e "$R1"
            R2=${dedup_Dir}/${ID}_deduped_R2.fastq.gz; echo -e "$R2"
            
            

        H2 "Output"
            out=${moduleDir}
            H3 "TrimGalore"
            echo ${out}/${ID}_trimmed_1.fastq.gz
            echo ${out}/${ID}_trimmed_2.fastq.gz

            H3 "FastQC"
            qc_out=${out}/QC_report
            if [[ ! -d ${qc_out} ]]; then mkdir -p ${qc_out}; fi
            echo ${qc_out}/${ID}_trimmed_1_fastqc.html
            echo ${qc_out}/${ID}_trimmed_1_fastqc.zip
            echo ${qc_out}/${ID}_trimmed_2_fastqc.html
            echo ${qc_out}/${ID}_trimmed_2_fastqc.zip

            if [[ ! -d ${out}/COMPLETE ]]; then mkdir -p $out/COMPLETE; fi

    H3 "[ Start ]"
    /bin/date
    SECONDS=0
    start=$SECONDS
    Complete_tag=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env

# Run functions ------------------------------------------------------------------------------------------------------------
    
    STEP="TrimGalore"
        H1 "$STEP"
        Trimmed_files=("${out}/${ID}_trimmed_1.fastq.gz" "${out}/${ID}_trimmed_2.fastq.gz")
        # echo ${Trimmed_files[@]}
        output_exists=$(test_for_output "${Trimmed_files[@]}")
        # echo $output_exists

        if ! $output_exists; then
            
            CMD="trim_galore --no_report_file --paired $R1 $R2 -o $out"; echo $CMD
            $CMD; if [[ $? -ne 0 ]]; then error "Something went wrong! Exiting..."; fi

            comment "Fixing default TrimGalore naming output"
            tmp=${R1%_*}; sample=${tmp##*/}
            mv ${out}/${sample}_R1_val_1.fq.gz ${out}/${ID}_trimmed_1.fastq.gz
            mv ${out}/${sample}_R2_val_2.fq.gz ${out}/${ID}_trimmed_2.fastq.gz

        fi
        step_completion "${Trimmed_files[@]}"


    STEP="QC Report on Trimmed Reads"

        H1 "$STEP"
        reads_1=${out}/${ID}_trimmed_1.fastq.gz
        reads_2=${out}/${ID}_trimmed_2.fastq.gz

        QC_files=("${qc_out}/${ID}_trimmed_1_fastqc.zip" "${qc_out}/${ID}_trimmed_1_fastqc.html")
        QC_files+=("${qc_out}/${ID}_trimmed_2_fastqc.zip" "${qc_out}/${ID}_trimmed_2_fastqc.html")
        echo ${QC_files[@]}
        output_exists=$(test_for_output "${QC_files[@]}")

        if ! $output_exists; then
            
            CMD="fastqc -q -t $SLURM_NTASKS -o ${qc_out} -f fastq $reads_1 $reads_2"; echo $CMD
            $CMD
            if [[ $? -ne 0 ]]; then error "Something went wrong! Exiting..."; fi

        fi
        step_completion "${QC_files[@]}"


    H1 "Completion"
        
        output_exists=$(test_for_output "${Complete_tag[@]}")
        if $output_exists; then touch ${out}/COMPLETE/$ID; fi

# Unload environments ------------------------------------------------------------------------------------------------------
    conda deactivate
    module unload anaconda3/2023.09

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
