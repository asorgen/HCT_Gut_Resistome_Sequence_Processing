#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Concatenates a sample's host-decontaminated paired reads into one FASTQ
    #    (HUMAnN has no native paired-end mode; concatenation is HUMAnN's own
    #    documented recommendation for paired-end input)
    # 2. Runs HUMAnN (nucleotide search vs ChocoPhlAn, translated search vs
    #    UniRef90, gene-family/pathway quantification) in one command
    #
    # Output: per-sample whole-community functional profile:
    #   {ID}_genefamilies.tsv   — UniRef90 gene-family abundance (RPK)
    #   {ID}_pathabundance.tsv  — MetaCyc pathway abundance
    #   {ID}_pathcoverage.tsv   — MetaCyc pathway coverage
    #
    # Companion analysis for [[Build taxonomy and gene-family companion
    # analyses for SNP paper]] — whole-community gene-family profiling to pair
    # with the null SNP-level in-lineage-selection result. NOTE: humann_opts
    # (Duke_short-read.config) is an unbenchmarked starting guess pending
    # single-sample validation — see that config's comment.
    #
    # Created on 2026-07-20
    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics
    # Version: 1

# Config files -------------------------------------------------------------------------------------------------------------
    source $pipelineConfig
    source $config_file
    source $bashrc
    source $bash_profile
    source $module_functions
    source $print_functions

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 5.8_humann.sh"
        echo -e "Whole-community gene-family/pathway profiling using HUMAnN v3.8."
        echo -e "1. Concatenate paired reads (R1+R2) into one FASTQ"
        echo -e "2. humann — nucleotide search (ChocoPhlAn) + translated search (UniRef90) + quantification"

    H1 "Job Context"
        print "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        print "Running on host: $(hostname)"
        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))
        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")
        print "----- Resources Requested -----"
        print "Nodes:           $SLURM_NNODES"
        print "Cores / node:    $SLURM_CPUS_PER_TASK"
        print "Total memory:    ${Total_Gb} GB"
        print "Wall-clock time: $JobTime"
        print "-------------------------------"

    H1 "Variables"
        comment "SampleID: ${ID}"

        H2 "Input"
            R1=${clean_readDir}/${ID}_1.fastq.gz
            R2=${clean_readDir}/${ID}_2.fastq.gz
            echo -e "$R1"
            echo -e "$R2"
            [[ ! -f "$R1" ]] && { echo "ERROR: Input file not found: $R1"; exit 1; }
            [[ ! -f "$R2" ]] && { echo "ERROR: Input file not found: $R2"; exit 1; }

        H2 "Output"
            if [[ ! -d ${moduleDir}/${ID} ]]; then mkdir -p ${moduleDir}/${ID}; fi
            concatFastq=${moduleDir}/${ID}/${ID}_concat.fastq.gz
            geneFamiliesOut=${moduleDir}/${ID}/${ID}_genefamilies.tsv
            pathAbundOut=${moduleDir}/${ID}/${ID}_pathabundance.tsv
            pathCovOut=${moduleDir}/${ID}/${ID}_pathcoverage.tsv
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=("$geneFamiliesOut" "$pathAbundOut")
    Intermediate_files=("$concatFastq")

# Load environments --------------------------------------------------------------------------------------------------------
    module load humann/3.8

# Run functions ------------------------------------------------------------------------------------------------------------

func="Concatenate paired reads"
    H1 "$func"
    if [[ -s "$concatFastq" ]]; then
        comment "Output already found. Skipping..."
    else
        start=$SECONDS
        cat "$R1" "$R2" > "$concatFastq"
        [[ $? -ne 0 ]] && error "Concatenation failed for $ID"
        end=$SECONDS; comment "Elapsed: $(elapsed_time $((end - start)))"
    fi
    step_completion "$concatFastq"


func="HUMAnN gene-family/pathway profiling"
    H1 "$func"
    if [[ -s "$geneFamiliesOut" ]]; then
        comment "Output already found. Skipping..."
    else
        start=$SECONDS
        humann \
            --input "$concatFastq" \
            --output "${moduleDir}/${ID}" \
            --output-basename "${ID}" \
            --threads $SLURM_CPUS_PER_TASK \
            --nucleotide-database "$HUMANN_CHOCOPHLAN" \
            --protein-database "$HUMANN_UNIREF" \
            --remove-temp-output
        [[ $? -ne 0 ]] && error "humann failed for $ID"
        N_GENEFAM=$(($(wc -l < "$geneFamiliesOut") - 1))
        comment "Gene families quantified: $N_GENEFAM"
        end=$SECONDS; comment "Elapsed: $(elapsed_time $((end - start)))"
    fi
    step_completion "$geneFamiliesOut" "$pathAbundOut"


# Completion ---------------------------------------------------------------------------------------------------------------
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then
        touch ${moduleDir}/COMPLETE/$ID
        for int_file in "${Intermediate_files[@]}"; do
            echo -e "rm $int_file"; rm -f "$int_file"
        done
    else
        error "Not all outputs were created."
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
