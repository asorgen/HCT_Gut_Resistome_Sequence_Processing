#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Unzips the final cleaned sequence files.
    # 2. Performs binning with MetaBat2 using metaWRAP's binning module.
    #     1. Indexes assembly file (bwa index).
    #     2. Aligns reads back to the assembly (bwa mem).
    #     3. Sorts the alignment file (samtools sort)
    #     4. Makes a contig depth file (jgi_summarize_bam_contig_depths).
    #     5. Performs binning with MetaBat2
    # 3. Performs binning with MaxBin2 using metaWRAP's binning module.
    #     1. Indexes assembly file (bwa index).
    #     2. Aligns reads back to the assembly (bwa mem).
    #     3. Sorts the alignment file (samtools sort).
    #     4. Makes a contig depth file (jgi_summarize_bam_contig_depths).
    #     5. Splits the contig depth file into multiple files.
    #     6. Performs binning with MaxBin2.
    # 4. Performs binning with CONCOCT using metaWRAP's binning module.
    #     1. Indexes assembly file (bwa index).
    #     2. Aligns reads back to the assembly (bwa mem).
    #     3. Sorts the alignment file (samtools sort).
    #     4. Indexes the bam alignment file.
    #     5. Cuts contigs into 10kb fragments for CONCOCT (cut_up_fasta.py).
    #     6. Estimates contig fragment coverage (concoct_coverage_table.py).
    #     7. Performs binning with CONCOCT.
    #     8. Merges 10kb fragments back into contigs (merge_cutup_clustering.py).
    #     9. Splits contigs into bins (split_concoct_bins.py).

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

# Set function for output comments -----------------------------------------------------------------------------------------
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Usage"
echo -e "This script does the following:"
echo -e "1. Unzips the final cleaned sequence files."
echo -e "2. Performs binning with MetaBat2 using metaWRAP's binning module."
echo -e "    1. Indexes assembly file (bwa index)."
echo -e "    2. Aligns reads back to the assembly (bwa mem)."
echo -e "    3. Sorts the alignment file (samtools sort)."
echo -e "    4. Makes a contig depth file (jgi_summarize_bam_contig_depths)."
echo -e "    5. Performs binning with MetaBat2"
echo -e "3. Performs binning with MaxBin2 using metaWRAP's binning module."
echo -e "    1. Indexes assembly file (bwa index)."
echo -e "    2. Aligns reads back to the assembly (bwa mem)."
echo -e "    3. Sorts the alignment file (samtools sort)."
echo -e "    4. Makes a contig depth file (jgi_summarize_bam_contig_depths)."
echo -e "    5. Splits the contig depth file into multiple files."
echo -e "    6. Performs binning with MaxBin2."
echo -e "4. Performs binning with CONCOCT using metaWRAP's binning module."
echo -e "    1. Indexes assembly file (bwa index)."
echo -e "    2. Aligns reads back to the assembly (bwa mem)."
echo -e "    3. Sorts the alignment file (samtools sort)."
echo -e "    4. Indexes the bam alignment file."
echo -e "    5. Cuts contigs into 10kb fragments for CONCOCT (cut_up_fasta.py)."
echo -e "    6. Estimates contig fragment coverage (concoct_coverage_table.py)."
echo -e "    7. Performs binning with CONCOCT."
echo -e "    8. Merges 10kb fragments back into contigs (merge_cutup_clustering.py)."
echo -e "    9. Splits contigs into bins (split_concoct_bins.py)."

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
        comment "SampleID (ID): ${ID}"

        H2 "Input"
            R1=${clean_readDir}/${ID}_1.fastq.gz; echo -e "${R1}"
            R2=${clean_readDir}/${ID}_2.fastq.gz; echo -e "${R2}"
            echo -e "${evaluationDir}/${ID}_final_assembly.fasta"

        H2 "Output"
            out=${moduleDir}/$ID
            if [[ ! -d ${out} ]]; then mkdir -p $out; fi

            metabat2_OUT=${moduleDir}/${ID}/metabat2_bins; echo -e "${metabat2_OUT}"
            maxbin2_OUT=${moduleDir}/${ID}/maxbin2_bins; echo -e "${maxbin2_OUT}"
            concoct_OUT=${moduleDir}/${ID}/concoct_bins; echo -e "${concoct_OUT}"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p $moduleDir/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env

STEP="Unzip the final cleaned sequence files."
    reads_1=${out}/${ID}_1.fastq
    reads_2=${out}/${ID}_2.fastq
    Unzip_files=(${reads_1} ${reads_2})
    output_exists=$(test_for_output "${Unzip_files[@]}")
    if ! "$output_exists"; then
        H1 "$STEP"
        gunzip -c $R1 > ${reads_1}
        gunzip -c $R2 > ${reads_2}
    fi 
    substep_completion "${Unzip_files[@]}"


H1 "Binning"
    H2 "MetaBat2"
        start=$SECONDS
        metawrap binning \
            -t $SLURM_CPUS_PER_TASK \
            -a ${evaluationDir}/${ID}_final_assembly.fasta \
            -o ${moduleDir}/${ID} \
            --metabat2 \
            ${reads_1} \
            ${reads_2}
        if [[ ! -d $metabat2_OUT ]]; then error "Something went wrong with MetaBat2. Exiting..."; fi
        end=$SECONDS; duration=$(( end-start ))
        comment "Time spent: $(elapsed_time "$duration")"

    H2 "MaxBin2"
        start=$SECONDS
        metawrap binning \
            -t $SLURM_CPUS_PER_TASK \
            -a ${evaluationDir}/${ID}_final_assembly.fasta \
            -o ${moduleDir}/${ID} \
            --maxbin2 \
            ${reads_1} \
            ${reads_2}
        if [[ ! -d $maxbin2_OUT ]]; then error "Something went wrong with MetaBat2. Exiting..."; fi
        end=$SECONDS; duration=$(( end-start ))
        comment "Time spent: $(elapsed_time "$duration")"

    H2 "CONCOCT"
        start=$SECONDS
        metawrap binning \
            -t $SLURM_CPUS_PER_TASK \
            -a ${evaluationDir}/${ID}_final_assembly.fasta \
            -o ${moduleDir}/${ID} \
            --concoct \
            ${reads_1} \
            ${reads_2}
        if [[ ! -d $metabat2_OUT ]]; then error "Something went wrong with MetaBat2. Exiting..."; fi
        end=$SECONDS; duration=$(( end-start ))
        comment "Time spent: $(elapsed_time "$duration")"

conda deactivate
module unload anaconda3/2023.09

# Completion status
    if [[ -d $metabat2_OUT && -d $maxbin2_OUT && -d $concoct_OUT ]]; then
        touch ${moduleDir}/COMPLETE/$ID

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm $int_file
        done

    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
