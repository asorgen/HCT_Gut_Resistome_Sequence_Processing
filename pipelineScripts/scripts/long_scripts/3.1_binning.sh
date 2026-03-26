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
    source $print_functions


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
        OMP_NUM_THREADS=$SLURM_NTASKS
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1000 ))

        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

        echo 
        print "----- Resources Requested -----"
        print "Nodes:            $SLURM_NNODES"
        print "Cores / node:     $SLURM_NTASKS"
        print "Total memory:     $Total_Gb Gb"
        print "Wall-clock time:  $JobTime"
        print "-------------------------------"

    H1 "Variables"
        comment -e "SampleID (ID): ${ID}"

        H2 "Input"
            long_reads=${clean_readDir}/${ID}_ont.fastq; echo -e $long_reads
            if [[ "$readType" == "hybrid" ]]; then
                    reads_1=${clean_readDir}/${ID}_1.fastq; echo -e $reads_1
                    reads_2=${clean_readDir}/${ID}_2.fastq; echo -e $reads_2
            fi
            ASSEMBLY=${evaluationDir}/${ID}_final_assembly.fasta; echo -e "${ASSEMBLY}"

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


# Default params & set up ------------------------------------------------------------------------------------------------
    threads=$SLURM_CPUS_PER_TASK; mem=$Total_Gb; len=1000; out=${moduleDir}/${ID}; s
    markers=107

    if [ $len -lt 1500 ]; then
        metabat_len=1500
    else
        metabat_len=$len
    fi

    if [[ ! -d ${out}/work_files ]]; then mkdir -p ${out}/work_files; fi
    cp $ASSEMBLY ${out}/work_files/assembly.fa


# Asssembly alignment ------------------------------------------------------------------------------------------------------
    STEP="Assembly alignment"
    
    bamFile=${out}/work_files/${ID}_long.bam
    output_exists=$(test_for_output "${bamFile}")
    if ! "$output_exists"; then
        H1 "$STEP"

        # module load minimap2
        # module load samtools
        # module load bwa

        # Index assembly
            STEP="Indexing assembly"
            outfile=${out}/work_files/filtered_assembly.mmi
            output_exists=$(test_for_output "${outfile}")
            if ! "$output_exists"; then
                H3 $STEP
                minimap2 -d ${out}/work_files/filtered_assembly.mmi ${out}/work_files/assembly.fa
                
                # if [[ "$readType" == "hybrid" ]]; then
                #     bwa index ${out}/work_files/assembly.fa
                #     if [[ $? -ne 0 ]]; then error "Something went wrong with bwa assembly indexing. Exiting"; fi
                # fi
            fi; substep_completion $outfile

        # Map contigs back to assembly
            STEP="Mapping contigs back to assembly"
            outfile=${out}/work_files/${ID}_long.sam
            output_exists=$(test_for_output "${outfile}")
            if ! "$output_exists"; then
                H3 $STEP
                minimap2 -a -x map-ont  ${out}/work_files/filtered_assembly.mmi $long_reads > ${out}/work_files/${ID}_long.sam
                
                # if [[ "$readType" == "hybrid" ]]; then
                #     bwa mem -v 1 -t $threads ${out}/work_files/assembly.fa $reads_1 $reads_2 > ${out}/work_files/${ID}_short.sam
                #     if [[ $? -ne 0 ]]; then error "Something went wrong with bwa read mapping. Exiting"; fi
                # fi
            fi; substep_completion ${outfile}


        # Sort the mapping file 
            STEP="Sorting the mapping file"
            outfile=${out}/work_files/${ID}_long_unsort.bam
            output_exists=$(test_for_output "${outfile}")
            if ! "$output_exists"; then
                samtools view -bS ${out}/work_files/${ID}_long.sam > ${out}/work_files/${ID}_long_unsort.bam
                samtools sort -@ $threads -o ${out}/work_files/${ID}_long.bam ${out}/work_files/${ID}_long_unsort.bam
                
                # if [[ "$readType" == "hybrid" ]]; then
                #     samtools view -bS ${out}/work_files/${ID}_short.sam > ${out}/work_files/${ID}_short_unsort.bam
                #     samtools sort -@ $threads -o ${out}/work_files/${ID}_short.bam ${out}/work_files/${ID}_short_unsort.bam
                #     if [[ $? -ne 0 ]]; then error "Something went wrong while sorting the short mapping file. Exiting"; fi
                # fi
            fi; substep_completion ${outfile}

        # module unload minimap2
        # module unload samtools
        # module unload bwa

    fi; substep_completion $bamFile
    exit 0

# MetaBat2 ---------------------------------------------------------------------------------------------------
    func="MetaBat2"
    H1 "Binning with $func"
    outputFile=${out}/metabat2_bins
    if [[ -d "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        H3 "Making contig depth file"
            jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/metabat_depth.txt ${out}/work_files/*.bam
            if [[ $? -ne 0 ]]; then error "Something went wrong while summarizing contig depths. Exiting"; fi

        H3 "Starting $func"
            # mkdir -p ${out}/metabat2_bins
            metabat2 -i $ASSEMBLY -a ${out}/work_files/metabat_depth.txt \
                 -o ${out}/metabat2_bins/bin -m $metabat_len -t $threads --unbinned
            if [[ $? -ne 0 ]]; then echo "Something went wrong while running $func. Exiting"; rm -r ${out}/metabat2_bins; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ $? -eq 0 ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi

# MaxBin2 ---------------------------------------------------------------------------------------------------
    func="MaxBin2"
    H1 "Binning with $func"
    outputFile=${out}/maxbin2_bins
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        H3 "Making contig depth file"
            jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/mb2_master_depth.txt --noIntraDepthVariance ${out}/work_files/*.bam
            if [[ $? -ne 0 ]]; then error "Something went wrong while summarizing contig depths. Exiting"; fi

        H3 "Splitting the master contig depth file into individual files for $func input"
            #calculate total numper of columns
            A=($(head -n 1 ${out}/work_files/mb2_master_depth.txt)) 
            N=${#A[*]}

            if [ -f ${out}/work_files/mb2_abund_list.txt ]; then rm ${out}/work_files/mb2_abund_list.txt; fi
            for i in $(seq 4 $N); do 
                sample=$(head -n 1 ${out}/work_files/mb2_master_depth.txt | cut -f $i)
                echo "Processing depth file..."
                grep -v totalAvgDepth ${out}/work_files/mb2_master_depth.txt | cut -f 1,$i > ${out}/work_files/mb2_${ID}.txt
                echo ${out}/work_files/mb2_${ID}.txt >> ${out}/work_files/mb2_abund_list.txt
            done
            if [[ $? -ne 0 ]]; then echo "Something went wrong while running $func. Exiting"; rm -r ${out}/metabat2_bins; fi


        H3 "Starting $func"
            mkdir ${out}/work_files/maxbin2_out
            run_MaxBin.pl -contig $ASSEMBLY -markerset $markers -thread $threads -min_contig_length $len \
                -out ${out}/work_files/maxbin2_out/bin \
                -abund_list ${out}/work_files/mb2_abund_list.txt

            if [[ $? -ne 0 ]]; then 
                echo "Something went wrong while running $func. Exiting"
            elif [[ $(ls ${out}/work_files/maxbin2_out/ | grep bin | grep .fasta | wc -l) -lt 1 ]]; then 
                echo "MaxBin2 did not pruduce a single bin. Something went wrong. Exiting."
            else
                mkdir ${out}/maxbin2_bins
                N=0
                for i in $(ls ${out}/work_files/maxbin2_out/ | grep bin | grep .fasta); do
                    cp ${out}/work_files/maxbin2_out/$i ${out}/maxbin2_bins/bin.${N}.fa
                    N=$((N + 1))
                done
                rm -r ${out}/work_files/maxbin2_out
            fi
            


        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ $? -eq 0 ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi

# CONCOCT ---------------------------------------------------------------------------------------------------
    func="CONCOCT"
    H1 "Binning with $func"
    outputFile=${out}/concoct_bins
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        if [[ ! -s ${out}/work_files/concoct_depth.txt ]]; then
            H3 "Indexing .bam alignment files..."
                samtools index -@ $threads -b ${out}/work_files/${ID}_long.bam
                if [[ $? -ne 0 ]]; then error "Something went wrong with indexing bam files. Exiting."; fi

            H3 "Cutting up contigs into 10kb fragments for CONCOCT..."
                cut_up_fasta.py $ASSEMBLY -c 10000 --merge_last -b ${out}/work_files/assembly_10K.bed -o 0 > ${out}/work_files/assembly_10K.fa
                if [[ $? -ne 0 ]]; then error "Something went wrong with cutting up contigs. Exiting."; fi

            H3 "Estimating contig fragment coverage..."   
                CMD="concoct_coverage_table.py ${out}/work_files/assembly_10K.bed ${out}/work_files/${ID}_long.bam > ${out}/work_files/concoct_depth.txt"
                $(eval $CMD)
                if [[ $? -ne 0 ]]; then error "Something went wrong with estimating fragment abundance. Exiting..."; fi
        fi


        H3 "Starting binning with CONCOCT..."
            mkdir -p ${out}/work_files/concoct_out
            # conda info --envs
            # conda activate concoct_env

            concoct -l $len -t $threads \
                --coverage_file ${out}/work_files/concoct_depth.txt \
                --composition_file ${out}/work_files/assembly_10K.fa \
                -b ${out}/work_files/concoct_out
            if [[ $? -ne 0 ]]; then error "Something went wrong with binning with CONCOCT. Exiting..."; fi

        H3 "Merging 10kb fragments back into contigs"
            merge_cutup_clustering.py ${out}/work_files/concoct_out/clustering_gt${len}.csv > ${out}/work_files/concoct_out/clustering_gt${len}_merged.csv
            if [[ $? -ne 0 ]]; then error "Something went wrong with merging fragments. Exiting..."; fi

        H3 "Splitting contigs into bins"
            mkdir ${out}/concoct_bins
            ${SOFT}/split_concoct_bins.py ${out}/work_files/concoct_out/clustering_gt${len}_merged.csv ${out}/work_files/assembly.fa ${out}/concoct_bins
            if [[ $? -ne 0 ]]; then echo "Something went wrong with splitting contigs into bins. Exiting..."; rm -r ${out}/concoct_bins; exit 1; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ $? -eq 0 ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
        conda deactivate

    fi

rm -r ${out}/work_files
conda deactivate

# Completion status
    metabat2_OUT=${moduleDir}/${ID}/metabat2_bins
    maxbin2_OUT=${moduleDir}/${ID}/maxbin2_bins
    concoct_OUT=${moduleDir}/${ID}/concoct_bins
    if [[ -d $metabat2_OUT && -d $maxbin2_OUT && -d $concoct_OUT ]]; thenm
        touch ${moduleDir}/${ID}/COMPLETE
    fi


H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
