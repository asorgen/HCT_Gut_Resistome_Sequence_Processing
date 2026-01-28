#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------

    # This script does the following:
    # 1. Unzips the final cleaned sequence files.
    # 2. Assembles the reads using metaWRAP's assembly module:
    #     1. Assembles reads using metaspades.
    #     2. Sorts out the unassembled reads.
    #         1. Removes scaffolds less than 1.5kb (rm_short_contigs.py).
    #         2. Indexes the remaining scaffolds (bwa index).
    #         3. Sorts out and stores reads that don't map back to the assembly (bwa mem, sam_to_fastq.py).
    #     3. Assembes the unassembled reads with megahit.
    #     4. Combines and formats the two assemblies.
    #     5. Performs a quality assessment of the final assembly with Quast.

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
    #-----------------------------------------------

    #SBATCH --mail-user=${email}

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
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }
    test_for_output() {
        local File_list=("$@")  # grab all args as one big array
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
    echo -e "1. Unzips the final cleaned sequence files."
    echo -e "2. Assembles the reads using metaWRAP's assembly module:"
    echo -e "    1. Assembles reads using metaspades."
    echo -e "    2. Sorts out the unassembled reads."
    echo -e "        1. Removes scaffolds less than 1.5kb (rm_short_contigs.py)."
    echo -e "        2. Indexes the remaining scaffolds (bwa index)."
    echo -e "        3. Sorts out and stores reads that don't map back to the assembly (bwa mem, sam_to_fastq.py)."
    echo -e "    3. Assembes the unassembled reads with megahit."
    echo -e "    4. Combines and formats the two assemblies."
    echo -e "    5. Performs a quality assessment of the final assembly with Quast."

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
        echo -e "${clean_readDir}/${ID}_1.fastq"
        echo -e "${clean_readDir}/${ID}_2.fastq"
        H2 "Output"
        echo -e "Assembly files will be deposited to ${moduleDir}/${ID}"
        if [[ ! -d ${moduleDir}/${ID} ]]; then mkdir ${moduleDir}/${ID}; fi
        if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi
    

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    start=$SECONDS

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    
# Run functions ------------------------------------------------------------------------------------------------------------
    Unzip_files=(${moduleDir}/${ID}/${ID}_1.fastq ${moduleDir}/${ID}/${ID}_2.fastq)
    Final_assembly=(${moduleDir}/${ID}/assembly.fasta)
    Final_assembly_exists=$(test_for_output "${Final_assembly[@]}")
    
    if ! $Final_assembly_exists; then
        STEP="Unzip cleaned sequence files"
            
            Unzip_files_exist=$(test_for_output "${Unzip_files[@]}")

            if ! $Unzip_files_exist; then
                H1 "$STEP"
                echo "gunzip -c ${clean_readDir}/${ID}_1.fastq.gz > ${moduleDir}/${ID}/${ID}_1.fastq"
                gunzip -c ${clean_readDir}/${ID}_1.fastq.gz > ${moduleDir}/${ID}/${ID}_1.fastq
                echo "gunzip -c ${clean_readDir}/${ID}_2.fastq.gz > ${moduleDir}/${ID}/${ID}_2.fastq"
                gunzip -c ${clean_readDir}/${ID}_2.fastq.gz > ${moduleDir}/${ID}/${ID}_2.fastq
            fi

            substep_completion "${Unzip_files[@]}"

        STEP="metaSPAdes Assembly"

            # if [[ -d ${moduleDir}/${ID} ]]; then rm -rf ${moduleDir}/${ID}; fi
            conda activate metawrap-env

            output_exists=$(test_for_output "${Final_assembly[@]}")
            if ! $output_exists; then
                H1 "$STEP"
                metawrap assembly \
                    -1 ${moduleDir}/${ID}/${ID}_1.fastq \
                    -2 ${moduleDir}/${ID}/${ID}_2.fastq \
                    -m $Total_Gb \
                    -t $SLURM_CPUS_PER_TASK \
                    --metaspades \
                    -o ${moduleDir}/${ID}
                if [[ $? -ne 0 ]]; then error "Something went wrong! Exiting..."; fi

                # Default metaWRAP assembly file: final_assembly.fasta
                mv ${moduleDir}/${ID}/final_assembly.fasta ${moduleDir}/${ID}/assembly.fasta

            fi

            # step_completion "${moduleDir}/${ID}/assembly.fasta"


            conda deactivate
            module unload anaconda3/2023.09
    fi
    step_completion "${Final_assembly[@]}"


    


module_completion
H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
