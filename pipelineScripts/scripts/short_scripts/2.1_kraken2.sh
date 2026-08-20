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
    source module_functions
    source print_functions


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
        comment "SampleID: ${ID}"

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

    # Bracken's k-mer distribution is built for a specific read length, so estimating
    # short reads against a distribution built for longer ones biases abundances.
    # Pick the longest available distribution that does not exceed this sample's mean
    # read length. Duke/UNC (150-151bp) resolve to database150mers and are therefore
    # byte-identical to previous runs; cohorts with shorter reads (e.g. PRJNA890666,
    # much of which is ~126bp) drop to the appropriate shorter distribution.
    DEFAULT_KMER_DISTRIB=${KRAKEN2_DB}/database150mers.kmer_distrib

    read_len=$(head -n 4000 "${reads_1}" 2>/dev/null \
               | awk 'NR%4==2 {s+=length($0); n++} END {if (n>0) printf "%d", s/n}')

    # Available distributions, ascending
    avail=()
    for f in ${KRAKEN2_DB}/database*mers.kmer_distrib; do
        [[ -e $f ]] || continue
        n=$(basename "$f"); n=${n#database}; n=${n%mers.kmer_distrib}
        [[ $n =~ ^[0-9]+$ ]] && avail+=("$n")
    done
    if [[ ${#avail[@]} -gt 0 ]]; then
        IFS=$'\n' avail=($(sort -n <<<"${avail[*]}")); unset IFS
    fi

    kmer_distrib=$DEFAULT_KMER_DISTRIB
    if [[ -n "$read_len" && "$read_len" -gt 0 && ${#avail[@]} -gt 0 ]]; then
        best=""
        for n in "${avail[@]}"; do
            [[ $n -le $read_len ]] && best=$n
        done
        # Reads shorter than every available distribution fall back to the smallest.
        [[ -z $best ]] && best=${avail[0]}
        [[ -s ${KRAKEN2_DB}/database${best}mers.kmer_distrib ]] && \
            kmer_distrib=${KRAKEN2_DB}/database${best}mers.kmer_distrib
    fi

    comment "Mean read length: ${read_len:-unknown} bp"
    comment "Bracken k-mer distribution (reads): $(basename ${kmer_distrib})"

    for file in ${k_out}/*.kreport; do
        # Assembly kreports are contigs rather than reads, so they stay on the
        # original distribution - this keeps full-mode Duke/UNC output unchanged.
        if [[ $(basename "$file") == *_assembly.kreport ]]; then
            kd=$DEFAULT_KMER_DISTRIB
        else
            kd=$kmer_distrib
        fi
        python3 $bracken/est_abundance.py -i $file -k ${kd} --level S -o ${file%.*}.bracken.out
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



