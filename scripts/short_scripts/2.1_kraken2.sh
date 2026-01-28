#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Performs taxonomic classifications using a modification of metaWRAP's kraken2 module:
    #     1. Unzips the final cleaned sequence files.
    #     2. Performs taxonomic classification of sequence reads with Kraken2.
    #     3. Performs taxonomic classification of the final assembly with Kraken2.
    #     4. Translates the resulting Kraken2 output files
    #     5. Generates kronagrams of the output files
    # 2. Estimates species-level abundances with Bracken

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
        echo -e "1. Performs taxonomic classifications using a modification of metaWRAP's kraken2 module:"
        echo -e "    1. Unzips the final cleaned sequence files."
        echo -e "    2. Performs taxonomic classification of sequence reads with Kraken2."
        echo -e "    3. Performs taxonomic classification of the final assembly with Kraken2."
        echo -e "    4. Translates the resulting Kraken2 output files"
        echo -e "    5. Generates kronagrams of the output files"
        echo -e "2. Estimates species-level abundances with Bracken"

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
        comment -e "SampleID (ID): ${ID}"

        H2 "Input"
            R1=${clean_readDir}/${ID}_1.fastq.gz; echo -e "${R1}"
            R2=${clean_readDir}/${ID}_2.fastq.gz; echo -e "${R2}"
            if [[ -s ${evaluationDir}/${ID}_final_assembly.fasta ]]; then ASM=${evaluationDir}/${ID}_final_assembly.fasta; echo -e "$ASM"; fi
            
        H2 "Output"
            k_out=${moduleDir}/${ID}
            if [[ ! -d ${k_out} ]]; then mkdir -p $k_out; fi
            b_out=$brackenDir
            echo -e "Kraken2 and Krona output will be deposited to ${k_out}/"
            echo -e "Bracken will be deposited to ${b_out}/"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p $moduleDir/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    
# Run functions ------------------------------------------------------------------------------------------------------------

STEP="Unzip the sequence files"
    reads_1=${k_out}/${ID}_1.fastq
    reads_2=${k_out}/${ID}_2.fastq
    Unzip_files=(${reads_1} ${reads_2})
    output_exists=$(test_for_output "${Unzip_files[@]}")
    if ! "$output_exists"; then
        H1 "$STEP"
        gunzip -c $R1 > ${reads_1}
        gunzip -c $R2 > ${reads_2}
    fi
    substep_completion "${Unzip_files[@]}"

STEP="Kraken2"
    
    H1 "$STEP"
    conda activate metawrap-env
    if [[ -s ${evaluationDir}/${ID}_final_assembly.fasta ]]; then
        metawrap kraken2_bracken -o ${k_out} -t $SLURM_CPUS_PER_TASK ${k_out}/${ID}_*fastq ${evaluationDir}/${ID}_final_assembly.fasta
    else
        metawrap kraken2_bracken -o ${k_out} -t $SLURM_CPUS_PER_TASK ${k_out}/${ID}_*fastq
    fi
    conda deactivate



STEP="Bracken"

    H1 "$STEP"
    for file in ${k_out}/*.kreport; do
        python3 $bracken/est_abundance.py -i $file -k ${KRAKEN2_DB}/database150mers.kmer_distrib --level S -o ${file%.*}.bracken.out
        if [[ $? -ne 0 ]]; then error "Something went wrong while running Bracken. Exiting..."; fi
    done

    
    if [[ -s ${evaluationDir}/${ID}_final_assembly.fasta ]]; then
        mkdir -p ${b_out}/assembly
        mv ${k_out}/*_assembly.bracken.out ${b_out}/assembly
        Complete_tag+=(${b_out}/assembly/${ID}_assembly.bracken.out)
    fi
    
    mkdir -p ${b_out}/sr
    mv ${k_out}/*.bracken.out ${b_out}/sr
    step_completion ${b_out}/sr/${ID}.bracken.out


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
    module unload anaconda3/2023.09

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"



