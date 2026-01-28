#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. If needed: Indexes the human genome reference database (bowtie2).
    # 2. Unzips the the deduplicated, trimmed sequence read files for processing.
    # 3. Aligns sample reads to human genome reference database (GRCh38) with bowtie2.
    # 4. Converts alignment from a sam to bam format (samtools view).
    # 5. Filters out reads that are mapped to the human genome (samtools view).
    # 6. Sorts the bam file so that paired reads are matched (samtools sort).
    # 7. Writes the paired reads to their respective fastq files (samtools fastq).
    # 8. Performs a quality assessment of the final deduplicated, trimmed, and decontaminated sequence reads using FastQC.

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
        echo -e "1. If needed: Indexes the human genome reference database (bowtie2)."
        echo -e "2. Unzips the the deduplicated, trimmed sequence read files for processing."
        echo -e "3. Aligns sample reads to human genome reference database (GRCh38) with bowtie2."
        echo -e "4. Converts alignment from a sam to bam format (samtools view)."
        echo -e "5. Filters out reads that are mapped to the human genome (samtools view)."
        echo -e "6. Sorts the bam file so that paired reads are matched (samtools sort)."
        echo -e "7. Writes the paired reads to their respective fastq files (samtools fastq)."
        echo -e "8. Performs a quality assessment of the final deduplicated, trimmed, and decontaminated sequence reads using FastQC."

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
            R1=${trimmed_Dir}/${ID}_trimmed_1.fastq.gz; echo -e "$R1"
            R2=${trimmed_Dir}/${ID}_trimmed_2.fastq.gz; echo -e "$R2"
            echo -e "${GRCh38_DB}"

        H2 "Output"
            out=${moduleDir}
            H3 "Host decontamination"
            out_reads_1=${out}/${ID}_1.fastq.gz; echo $out_reads_1
            out_reads_2=${out}/${ID}_2.fastq.gz; echo $out_reads_2
            
            H3 "FastQC"
            qc_out=${out}/QC_report
            if [[ ! -d ${qc_out} ]]; then mkdir -p ${qc_out}; fi
            echo ${qc_out}/${ID}_cleaned_R1_fastqc.html
            echo ${qc_out}/${ID}_cleaned_R1_fastqc.zip
            echo ${qc_out}/${ID}_cleaned_R2_fastqc.html
            echo ${qc_out}/${ID}_cleaned_R2_fastqc.zip

            if [[ ! -d ${out}/COMPLETE ]]; then mkdir -p $out/COMPLETE; fi

    H3 "[ Start ]"
    /bin/date
    SECONDS=0
    start=$SECONDS
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env
    module load bowtie2
    module load samtools

# Run functions ------------------------------------------------------------------------------------------------------------
    

    STEP="Index human database"
        if [[ ! -d ${GRCh38_DB}/bowtie2_index ]]; then
            H1 "$STEP"
            mkdir ${GRCh38_DB}/bowtie2_index
            bowtie2-build ${GRCh38_DB}/GRCh38.p14.fa ${GRCh38_DB}/bowtie2_index/GRCh38.p14
            if [[ $? -ne 0 ]]; then error "Something went wrong! Exiting..."; fi
        fi

    Final_files=("${out}/${ID}_1.fastq.gz" "${out}/${ID}_2.fastq.gz")
    QC_files=("${qc_out}/${ID}_1_fastqc.zip" "${qc_out}/${ID}_1_fastqc.html")
    QC_files+=("${qc_out}/${ID}_2_fastqc.zip" "${qc_out}/${ID}_2_fastqc.html")

    all_output_exists=$(test_for_output "${Final_files[@]}" "${QC_files[@]}")
    
    if ! $all_output_exists; then

        final_fastq_files_exist=$(test_for_output "${Final_files[@]}")
        if ! $final_fastq_files_exist; then
            STEP="Unzip and rename the sequence files"
                reads_1=${out}/${ID}_trimmed_1.fastq
                reads_2=${out}/${ID}_trimmed_2.fastq
                Unzip_files=(${reads_1} ${reads_2})
                output_exists=$(test_for_output "${Unzip_files[@]}")
                if ! "$output_exists"; then
                    H1 "$STEP"
                    gunzip -c $R1 > ${reads_1}
                    gunzip -c $R2 > ${reads_2}
                fi
                substep_completion "${Unzip_files[@]}"
           
            STEP="Align reads to human database"
                sam_file=(${out}/${ID}_all.sam)
                output_exists=$(test_for_output "${sam_file[@]}")
                if ! $output_exists; then
                    H1 "$STEP"
                    bowtie2 -p $SLURM_NTASKS -x ${GRCh38_DB}/bowtie2_index/GRCh38.p14 -1 $reads_1 -2 $reads_2 -S ${out}/${ID}_all.sam
                fi
                substep_completion "${sam_file[@]}"

            STEP="Convert .sam file to .bam format"
                bam_file=(${out}/${ID}_all.bam)
                output_exists=$(test_for_output "${bam_file[@]}")
                if ! $output_exists; then
                    H1 "$STEP"
                    samtools view -bS ${out}/${ID}_all.sam > ${out}/${ID}_all.bam
                fi
                substep_completion "${bam_file[@]}"

            STEP="Filter for only unmapped reads"
                unmapped_bam_file=(${out}/${ID}_unmapped.bam)
                output_exists=$(test_for_output "${unmapped_bam_file[@]}")
                if ! $output_exists; then
                    H1 "$STEP"
                    samtools view -b -f 12 -F 256 ${out}/${ID}_all.bam > ${out}/${ID}_unmapped.bam
                fi
                substep_completion "${unmapped_bam_file[@]}"

            STEP="Sort .bam file so that paired reads are matched"
                unmapped_sorted_bam_file=(${out}/${ID}_unmapped_sorted.bam)
                output_exists=$(test_for_output "${unmapped_sorted_bam_file[@]}")
                if ! $output_exists; then
                    H1 "$STEP"
                    Gb_per_thread=$(( Total_Gb / SLURM_NTASKS ))
                    samtools sort -n -m ${Gb_per_thread}G -@ $SLURM_NTASKS ${out}/${ID}_unmapped.bam -o ${out}/${ID}_unmapped_sorted.bam
                fi
                substep_completion "${unmapped_sorted_bam_file[@]}"

            STEP="Write reads to respective .fastq.gz files"
                fastq_files_to_rename=(${out}/${ID}_cleaned_R1.fastq.gz ${out}/${ID}_cleaned_R2.fastq.gz)
                output_exists=$(test_for_output "${fastq_files_to_rename[@]}")
                if ! $output_exists; then
                    H1 "$STEP"
                    samtools fastq -@ $SLURM_NTASKS ${out}/${ID}_unmapped_sorted.bam -1 ${out}/${ID}_cleaned_R1.fastq.gz -2 ${out}/${ID}_cleaned_R2.fastq.gz -0 /dev/null -s /dev/null -n
                fi
                substep_completion "${fastq_files_to_rename[@]}"

            H1 "Copy and rename final fastq files"
                CMD="cp ${out}/${ID}_cleaned_R1.fastq.gz ${out}/${ID}_1.fastq.gz"; echo $CMD; $CMD
                CMD="cp ${out}/${ID}_cleaned_R2.fastq.gz ${out}/${ID}_2.fastq.gz"; echo $CMD; $CMD
            
        fi


        STEP="Post-QC Report"


            output_exists=$(test_for_output "${QC_files[@]}")
            
            reads_1=${out}/${ID}_1.fastq.gz
            reads_2=${out}/${ID}_2.fastq.gz

            if ! $output_exists; then
                
                H1 "$STEP"
                CMD="fastqc -q -t $SLURM_NTASKS -o ${qc_out} -f fastq $reads_1 $reads_2"; echo $CMD
                $CMD

            fi


    else
        step_completion "${Final_files[@]}"
        step_completion "${QC_files[@]}"
    fi




H1 "Completion"
    
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then 
        touch ${out}/COMPLETE/$ID

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            rm $int_file
        done
    fi

    # rm ${raw_readDir}/${ID}_1.fastq ${raw_readDir}/${ID}_2.fastq
    # rm ${dedup_Dir}/${ID}_deduped_R1.fastq.gz ${dedup_Dir}/${ID}_deduped_R1.fastq.gz
    # rm ${trimmed_Dir}/${ID}_trimmed_1.fastq ${trimmed_Dir}/${ID}_trimmed_2.fastq


# Unload environments ------------------------------------------------------------------------------------------------------
    conda deactivate
    module unload anaconda3/2023.09
    module unload bowtie2
    module unload samtools

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
