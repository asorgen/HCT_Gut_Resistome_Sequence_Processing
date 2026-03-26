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
    source $module_functions
    source $print_functions


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
        print "----- Resources Requested -----"
        print "Nodes:            $SLURM_NNODES"
        print "Cores / node:     $SLURM_CPUS_PER_TASK"
        print "Total memory:     $Total_Gb Gb"
        print "Wall-clock time:  $JobTime"
        print "-------------------------------"

    H1 "Variables"
        echo -e "SampleID (ID): ${ID}"
        H2 "Input"
            ONT=${clean_readDir}/${ID}_ont.fastq
            R1=${clean_readDir}/${ID}_1.fastq
            R2=${clean_readDir}/${ID}_2.fastq
            ASM=${evaluationDir}/${ID}_final_assembly.fasta

            if [[ "$readType" == "long_read" ]]; then
                echo -e "${ONT}"
                echo -e "$ASM"
            fi

            if [[ "$readType" == "hybrid" ]]; then
                echo -e "$R1"
                echo -e "$R2"
                echo -e "$ONT"
                echo -e "$ASM"
            fi

        H2 "Output"
            k_out=${moduleDir}/${ID}; if [[ ! -d ${k_out} ]]; then mkdir -p $k_out; fi
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

    H1 "Kraken2"
        start=$SECONDS

        conda activate metawrap-env

        if [[ "$readType" == "hybrid" ]]; then
            sr_out=${k_out}/${ID}.krak2
            sr_rep=${k_out}/${ID}.kreport

            
            if [ ! -s $sr_out ] || [ ! -s $sr_rep ]; then
                H2 "Processing short reads"
                output=$sr_out; report=$sr_rep
                CMD="kraken2 --use-names --db ${KRAKEN2_DB} --paired --threads $SLURM_CPUS_PER_TASK --output $output --report $report $R1 $R2"
                echo $CMD
                $CMD
            fi

            if [[ $? -ne 0 ]] || [[ ! -s $sr_out ]] || [[ ! -s $sr_rep ]]; then error "Something went wrong with Kraken2 while processing the short-read fastqs. Exiting..."; fi
        fi

        ont_out=${k_out}/${ID}_ONT.krak2
        ont_rep=${k_out}/${ID}_ONT.kreport

        asm_out=${k_out}/${ID}_assembly.krak2
        asm_rep=${k_out}/${ID}_assembly.kreport

        
        if [ ! -s $ont_out ] || [ ! -s $ont_rep ]; then
            H2 "Processing ONT reads"
            output=$ont_out; report=$ont_rep; input=$ONT
            CMD="kraken2 --use-names --db ${KRAKEN2_DB} --threads $SLURM_CPUS_PER_TASK --output $output --report $report $input"
            echo $CMD
            $CMD
            if [[ $? -ne 0 ]] || [[ ! -s $ont_out ]] || [[ ! -s $ont_rep ]]; then error "Something went wrong with Kraken2 while processing the ONT fastq. Exiting..."; fi
        fi


        if [ ! -s $asm_out ] || [ ! -s $asm_rep ]; then
            H2 "Processing assembly"
            output=$asm_out; report=$asm_rep; input=$ASM
            CMD="kraken2 --use-names --db ${KRAKEN2_DB} --threads $SLURM_CPUS_PER_TASK --output $output --report $report $input"
            echo $CMD
            $CMD
            if [[ $? -ne 0 ]] || [[ ! -s $asm_out ]] || [[ ! -s $asm_rep ]]; then error "Something went wrong with Kraken2 while processing the assembly fasta. Exiting..."; fi
        fi

        
    H1 "Kraken2 Translate"
        for file in ${k_out}/*.krak2; do
            comment "Translating $file"
            ${SOFT}/kraken2_translate.py ${KRAKEN2_DB} $file ${file%.*}.kraken2
            if [[ $? -ne 0 ]]; then error "Something went wrong with kraken2-translate. Exiting..."; fi
        done

    H1 "Kronogram"
        for file in ${k_out}/*.kraken2; do
            ${SOFT}/kraken_to_krona.py $file > ${file%.*}.krona
            if [[ $? -ne 0 ]]; then error "Something went wrong while making the krona file. Exiting..."; fi
        done

        ktImportText -o ${k_out}/kronagram.html ${k_out}/*krona
        if [[ ! -s ${k_out}/kronagram.html ]]; then error "Something went wrong while running KronaTools to make kronagram. Exiting..."; fi

        conda deactivate


    H1 "Bracken"

        for file in ${k_out}/*.kreport; do
            python3 ${HOME}/PROGRAMS/Bracken-2.7/src/est_abundance.py -i $file -k ${KRAKEN2_DB}/database150mers.kmer_distrib --level S -o ${file%.*}.bracken.out
            if [[ $? -ne 0 ]]; then error "Something went wrong while running Bracken. Exiting..."; fi
        done

        mkdir -p ${b_out}/assembly
        mv ${k_out}/*_assembly.bracken.out ${b_out}/assembly
        step_completion ${b_out}/assembly/${ID}_assembly.bracken.out
        
        mkdir -p ${b_out}/ont
        mv ${k_out}/*_ONT.bracken.out ${b_out}/ont
        step_completion ${b_out}/ont/${ID}_ONT.bracken.out

        if [[ "$readType" == "hybrid" ]]; then
            mkdir -p ${b_out}/sr
            mv ${k_out}/*.bracken.out ${b_out}/sr
            step_completion ${b_out}/sr/${ID}.bracken.out
        fi


        module unload anaconda3/2023.09

# Completion status
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then 
        touch ${moduleDir}/COMPLETE/$ID

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm $int_file
        done
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"



