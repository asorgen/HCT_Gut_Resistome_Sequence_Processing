#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script verifies all pipeline outputs exist and removes intermediate files for finished ONT pipeline samples.

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


source ${HOME}/.bashrc
source $config_file
source $pipelineConfig

# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }
    test_for_output() {
        File_list=("$@")
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
        File_list=("$@")

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

H1 "Usage"
    comment "This script verifies all pipeline outputs exist and removes intermediate files."

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

    H2 "Output"
        out=${datasetDir}/COMPLETE
        mkdir -p ${out}
        echo ${out}/${ID}


H3 "[ Start ]"
/bin/date
SECONDS=0
start=$SECONDS
Complete_tag=()
Intermediate_files=()

H1 "Double check outputs exist"

# 0.0_raw_reads — raw ONT file (intermediate, to be removed)
Intermediate_files+=(${raw_readDir}/${ID}.fastq.gz)

# 0.1_pre_qc

# 0.3_sequence_trim
Intermediate_files+=(${trimmed_Dir}/${ID}_trimmed.fastq.gz)

# 0.4_host_decontamination
STEP="Sequence Processing (ONT)"
Final_Fastq=(${clean_readDir}/${ID}_ont.fastq)
step_completion "${Final_Fastq[@]}"

# 2.1_kraken2 / 2.2_bracken
STEP="Taxonomic Classification"
step_completion "${brackenDir}/ont/${ID}_ONT.bracken.out"

# 5.4_rgi_bwt
if $run_rgi_bwt; then
    STEP="RGI BWT"
    step_completion "${rgi_bwt_dir}/kma_output/${ID}.rgi_kma.txt"
fi

# 1.2_evaluation (full mode)
if $run_eval; then
    STEP="Assembly Evaluation"
    step_completion "${evaluationDir}/${ID}_final_assembly.fasta"
    Intermediate_files+=(${evaluationDir}/${ID}_draft_assembly.fasta)
fi

# 3.3_reassemble_bins (full mode)
if $run_reassem; then
    STEP="Bin Reassembly"
    step_completion "${reassemDir}/${ID}/reassembled_bins.stats"
fi

# 5.5_AA_amr_assembly (full mode)
if $run_amr_aa_asm; then
    STEP="AA AMR from Assembly"
    step_completion "${datasetDir}/5.5_AA_amr_assembly/${ID}/${ID}_gene_annotations.tsv"
fi

H1 "Completion"

    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then
        CMD="touch ${out}/${ID}"; echo $CMD
        $CMD

        # Remove intermediate files (single files)
        for int_file in "${Intermediate_files[@]}"; do
            if [[ -e $int_file ]]; then
                CMD="rm $int_file"; echo $CMD
                $CMD
            fi
        done

        # Remove large intermediate directories (full mode only)
        if $full; then
            for int_dir in \
                "${assemblyDir}/${ID}" \
                "${binningDir}/${ID}" \
                "${refinedbinDir}/${ID}/work_files"; do
                if [[ -d $int_dir ]]; then
                    CMD="rm -rf $int_dir"; echo $CMD
                    $CMD
                fi
            done
        fi
    else
        error "Pipeline failed"
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
