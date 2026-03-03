#!/bin/bash

# Description
    #-----------------------------------------------
    # This script is used to bin the co-assembly using three different algorithms with the metaWRAP binning module.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 7, 2024

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1
    #-----------------------------------------------

# Slurm Resource Options

    #-----------------------------------------------
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

##SBATCH --mail-user=${email}

source ${HOME}/.bashrc
config_file=$(which config-metawrap)
source $config_file


# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }

H1 "Usage"
    comment "This script is used to bin the co-assembly using three different algorithms (MetaBat2, MaxBin2, and CONCOCT) with the metaWRAP binning module."

H1 "Job Context"
    OMP_NUM_THREADS=$SLURM_NTASKS
    comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
    comment "Running on host: `hostname`"

    Total_Gb=$(( SLURM_MEM_PER_NODE / 1000 ))

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
    echo -e "${evaluationDir}/${ID}_filtered_assembly.fasta"
    H2 "Output"
    echo -e "Binning output will be deposited to ${moduleDir}/${ID} under metabat2_bins/, maxbin2_bins/, and concoct_bins/."


H2 "[ Start ]"
/bin/date
SECONDS=0

module load anaconda3/2023.09
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
conda init
conda activate metawrap-env

# default params
threads=$SLURM_NTASKS; mem=$Total_Gb; len=1000; out=${moduleDir}/${ID}; ASSEMBLY=${evaluationDir}/${ID}_final_assembly.fasta
markers=107

if [ $len -lt 1500 ]; then
    metabat_len=1500
else
    metabat_len=$len
fi

mkdir -p ${out}/work_files
cp $ASSEMBLY ${out}/work_files/assembly.fa
long_reads=${clean_readDir}/${ID}_ont.fastq
reads_1=${clean_readDir}/${ID}_1.fastq
reads_2=${clean_readDir}/${ID}_2.fastq

func="Assembly alignment"
    H1 "$func"
    outputFile=${out}/work_files/${ID}_long.bam
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        # module load minimap2
        # module load samtools
        # module load bwa
        

        H3 "Indexing assembly"
            
            if [[ ! -s "${out}/work_files/filtered_assembly.mmi" ]]; then
                minimap2 -d ${out}/work_files/filtered_assembly.mmi ${out}/work_files/assembly.fa
                if [[ $? -ne 0 ]]; then error "Something went wrong with minimap2 assembly indexing. Exiting"; fi
                
                if [[ "$readType" == "hybrid" ]]; then
                    bwa index ${out}/work_files/assembly.fa
                    if [[ $? -ne 0 ]]; then error "Something went wrong with bwa assembly indexing. Exiting"; fi
                fi
            fi

        H3 "Mapping contigs back to assembly"
            
            if [[ ! -s "${out}/work_files/${ID}_long.sam" ]]; then
                minimap2 -a -x map-ont  ${out}/work_files/filtered_assembly.mmi $long_reads > ${out}/work_files/${ID}_long.sam
                if [[ $? -ne 0 ]]; then error "Something went wrong with minimap2 read mapping. Exiting"; fi
                
                if [[ "$readType" == "hybrid" ]]; then
                    bwa mem -v 1 -t $threads ${out}/work_files/assembly.fa $reads_1 $reads_2 > ${out}/work_files/${ID}_short.sam
                    if [[ $? -ne 0 ]]; then error "Something went wrong with bwa read mapping. Exiting"; fi
                fi
            fi


        H3 "Sorting the mapping file"
            
            if [[ ! -s "${out}/work_files/${ID}_long.bam" ]]; then
                samtools view -bS ${out}/work_files/${ID}_long.sam > ${out}/work_files/${ID}_long_unsort.bam
                samtools sort -@ ${SLURM_NTASKS} -o ${out}/work_files/${ID}_long.bam ${out}/work_files/${ID}_long_unsort.bam
                if [[ $? -ne 0 ]]; then error "Something went wrong while sorting the long mapping file. Exiting"; fi
                
                if [[ "$readType" == "hybrid" ]]; then
                    samtools view -bS ${out}/work_files/${ID}_short.sam > ${out}/work_files/${ID}_short_unsort.bam
                    samtools sort -@ ${SLURM_NTASKS} -o ${out}/work_files/${ID}_short.bam ${out}/work_files/${ID}_short_unsort.bam
                    if [[ $? -ne 0 ]]; then error "Something went wrong while sorting the short mapping file. Exiting"; fi
                fi
            fi

        # module unload minimap2
        # module unload samtools
        # module unload bwa


        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ $? -eq 0 ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi
    rm ${out}/work_files/*_unsort.bam

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
    if [[ -d $metabat2_OUT && -d $maxbin2_OUT && -d $concoct_OUT ]]; then
        touch ${moduleDir}/${ID}/COMPLETE
    fi


H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
