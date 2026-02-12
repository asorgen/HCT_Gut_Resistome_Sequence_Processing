#!/bin/bash

# Description
    #-------------------------------------------------------------------------------------------------------------------
    # This script is used to begin the HCT ARG Assembly pipeline for the Oxford Nanopore/Illumina (hybrid) assemblies.

    # It is important to note that it has been designed for a specific working directory. Therefore, the reproduction of the results will require small modifications of the script or the adaptation of your working directory.

    # Created on Nov 6, 2024

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 3

    # Required tools:
        # 1. metaWRAP metagenomic wrapper suite (https://github.com/bxlab/metaWRAP)
        # 2. Porechop ONT adapter trimmer (https://github.com/rrwick/Porechop)
        # 3. Flye de novo PacBio/ONT genome assembler (https://github.com/mikolmogorov/Flye)
        # 4. BBMap bioinformatic tools (https://github.com/BioInfoTools/BBMap)
        # 5. AMRFinder+ resistance gene identifier (https://github.com/ncbi/amr)
        # 6. RGI resistance gene identifier (https://github.com/arpcard/rgi)


    # This pipeline requires metaWRAP or a metawrap environtment for modules 00, 03-08. 

    # Ensure that the metaWRAP config file include paths to the necessary databases needed for this pipeline.
    # Instructions found here: https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md

    # Database      Size     Used in module
    # -----------------------------------------
    # CheckM        1.4Gb    binning, bin_refinement, reassemble_bins
    # Kraken2       125Gb    kraken2
    # NCBI_nt       71Gb     classify_bins
    # NCBI_tax      283Mb    classify_bins
    # Indexed hg38  20Gb     read_qc

    # Module 00 also requires Porechop.

    # Module 01 requires Flye.

    # Module 02 requires bbmap

    # Module 09 requires the tools AMRFinder+ (https://github.com/ncbi/amr) and RGI (https://github.com/arpcard/rgi) with the required databases.

    # The sampleList is a text file of the sample names in the following format:
    # #SampleID
    # D20248PRE
    # D20248D1

    # This pipeline was originally run on Red Hat Enterprise Linux 9.2 (Plow) using the Slurm Workload Manager.

    #-------------------------------------------------------------------------------------------------------------------

cohort=Duke
read=hybrid
dataset=${cohort}_${read}
export pipelineConfig=${dataset}-read.config
source pipelineScripts/configs/${pipelineConfig}
source $bashrc

if [[ ! -d $datasetDir ]]; then mkdir -p $datasetDir; fi
cd $datasetDir
mkdir -p LOGs

export long_seqPath
export short_seqPath
export readType


# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }
    job_lookup() { squeue -u asorgen --format='%.18i   %.9P   %.40j   %.1T   %.12M   %.10l   %.6D   %R' | awk -v id="$jobID" 'match($3,id) {print $1}'; }
    module_setup() {
        script=$1
        
        moduleDir="${script%.sh}"

        if [ ! -z "$2" ]; then
            moduleDir=$moduleDir/$2

        fi

        if [[ ! -d ${moduleDir} ]]; then mkdir -p $moduleDir; fi
        
        logDir=${moduleDir}/logs
        if [[ ! -d ${logDir} ]]; then mkdir -p $logDir; fi

        Complete_tag=${moduleDir}/${ID}/COMPLETE

        jobID="${script%%_*}"
        jobID=${jobID}_${ID}_${dataset}    
    }
    first_ID() {
        if [[ $module == 0 ]]; then 
            H2 "${ID}"
            ((module++))
        fi    
    }
    run_module() {
        
        # If the completion file exists
        if [[ -a "$Complete_tag" ]]; then
            
            Current_Job=COMPLETE

        else

            queued=$(job_lookup $jobID)
            # If the job is not queued
            if [ -z "$queued" ]; then

                # If the dependent job is complete
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then

                    first_ID
                    H3 "$header3"
                    Current_Job=$(sbatch \
                        $hpc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $Current_Job
                    echo "[ Log file ] -> ${dataset}/${logDir}/${ID}.${Current_Job##* }.log"
                    ((count++))

                else
                    H3 "$header3"
                    Current_Job=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $hpc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $Current_Job
                    echo "[ Log file ] -> ${dataset}/${logDir}/${ID}.${Current_Job##* }.log"
                    ((count++))

               fi # if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]

            else
                Current_Job=$queued
           fi # if [ -z "$queued" ]

        fi     
    }



if [[ ! -s ${dataset}_pipeline_${version}.out ]]; then
    H1 "Usage"
        comment "export version=$version"
        comment "nohup sh ./pipelineScripts/${dataset}_pipeline.sh >> ${dataset}/LOGs/${dataset}_pipeline_${version}.out 2>&1 &"
        
    H1 "Variables"
        comment "[ Raw short-read sequence directory ]: ${short_seqPath}"
        comment "[ Raw long-read sequence directory ]: ${long_seqPath}"

        
fi


##- Copy ONT raw sequences
    # Count the number of .fastq.gz files in the directory
    file_count=$(find "../Duke_long/raw_reads/" -maxdepth 1 -type f -name "*.fastq.gz" | wc -l)

    # Check if the count is 14
    if [ "$file_count" -eq 14 ]; then
        RAW_READS_JOB=COMPLETE
    else
        H1 "Copy raw nanopore files"
        comment "Copying raw Nanopore sequenece files."
        mkdir -p ../Duke_long/raw_reads
        script=cp_raw_ONT.sh
        RAW_READS_JOB=$(sbatch \
                --partition=Orion --nodes=1 --ntasks-per-node=16 --mem-per-cpu=2GB --time=24:00:00 \
                --job-name=cp_ONT -o ../Duke_long/LOGs/cp_raw_ONT.%A.log ${ROOT}/pipelineScripts/scripts/long_scripts/helper_scripts/${script})
        echo $RAW_READS_JOB
        echo "[ Log file ] -> ../Duke_long/LOGs/cp_raw_ONT.${RAW_READS_JOB##* }.log"
    fi


H1 "Submit Jobs"
# Initialize count variable
count=0

export totalSamples=`tail -n +2 $sampleList | wc -l`

for ID in $(tail -n +2 $sampleList); do
    
    module=0
    export ID

    ##- Read QC - short reads
        script=read_qc.sh
        moduleDir="${script%.sh}"
        export clean_readDir=clean_reads

        # If *_1.fastq and *_2.fastq are in HYBRID/clean_reads
        if [ -s "${clean_readDir}/${ID}_1.fastq" ] && [ -s "${clean_readDir}/${ID}_2.fastq" ]; then
            SHORT_QC_JOB=COMPLETE
        # If *_1.fastq and *_2.fastq are NOT in HYBRID/clean_reads
        else
            mkdir -p ${clean_readDir}

            short_clean_OUT1=../Duke_short/${clean_readDir}/${ID}_1.fastq
            short_clean_OUT2=../Duke_short/${clean_readDir}/${ID}_2.fastq

            qc_OUT1=../Duke_short/${moduleDir}/post-QC_reports/${ID}_1_fastqc.html
            qc_OUT2=../Duke_short/${moduleDir}/post-QC_reports/${ID}_2_fastqc.html

            # If *_1.fastq and *_2.fastq are in Duke_short/clean_reads AND *_1_fastqc.html and *_2_fastqc.html are in Duke_short/read_qc/post-QC_reports
            if [ -s "$short_clean_OUT1" ] && [ -s "$short_clean_OUT2" ] && [ -s "$qc_OUT1" ] && [ -s "$qc_OUT2" ]; then
                ln ${DATA_ROOT}/Duke_short/${clean_readDir}/${ID}_1.fastq ${DATA_ROOT}/Duke_hybrid/${clean_readDir}/${ID}_1.fastq
                ln ${DATA_ROOT}/Duke_short/${clean_readDir}/${ID}_2.fastq ${DATA_ROOT}/Duke_hybrid/${clean_readDir}/${ID}_2.fastq
                SHORT_QC_JOB=COMPLETE

            # If *_1.fastq and *_2.fastq are NOT in Duke_short/clean_reads AND *_1_fastqc.html and *_2_fastqc.html are NOT in Duke_short/read_qc/post-QC_reports
            else

                jobID=$(printf "%02d" $module)
                jobID=${jobID}_${ID}_short

                export moduleDir=../Duke_short/${moduleDir}
                logDir=${moduleDir}/${ID}
                mkdir -p $logDir
                mkdir -p ../Duke_short/${clean_readDir}
                mkdir -p ../Duke_short/${moduleDir}/pre-QC_reports
                mkdir -p ../Duke_short/${moduleDir}/post-QC_reports
                export R1_ext
                export R2_ext

                queued=$(job_lookup $jobID)
                # If this job is not already in the queue
                if [ -z "$queued" ]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Short-read QC"
                    SHORT_QC_JOB=$(sbatch \
                        $short_read_qc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${ROOT}/pipelineScripts/scripts/short_scripts/${script})
                    echo $SHORT_QC_JOB
                    echo "[ Log file ] -> ${logDir}/${ID}.${SHORT_QC_JOB##* }.log"
                # If this job is already in the queue
                else
                    SHORT_QC_JOB=$queued
                fi # If this job is not already in the queue

            fi # If *_1.fastq and *_2.fastq are in Duke_short/clean_reads

        fi # If *_1.fastq and *_2.fastq are in HYBRID/clean_reads

    ##- Read QC - long reads

        script=read_qc.sh
        moduleDir="${script%.sh}"
        export clean_readDir=clean_reads

        DEPENDENT_JOB=$RAW_READS_JOB

        # If *_ont.fastq is already in HYBRID/clean_reads
        if [[ -s "${clean_readDir}/${ID}_ont.fastq" ]]; then
            LONG_QC_JOB=COMPLETE

        # If *_ont.fastq is NOT already in HYBRID/clean_reads
        else
            long_clean_OUT=../Duke_long/${clean_readDir}/${ID}_ont.fastq
            qc_OUT=../Duke_long/${moduleDir}/post-QC_reports/${ID}_fastqc.html

            # If *_ont.fastq is already in Duke_long/clean_reads AND *_fastqc.html is in Duke_long/read_qc/post-QC_reports
            if [ -s "$long_clean_OUT" ] && [ -s "$qc_OUT" ]; then
                ln ${DATA_ROOT}/Duke_long/${clean_readDir}/${ID}_ont.fastq ${DATA_ROOT}/Duke_hybrid/${clean_readDir}/${ID}_ont.fastq
                LONG_QC_JOB=COMPLETE

            # If *_ont.fastq is NOT already in Duke_long/clean_reads AND *_fastqc.html is NOT in Duke_long/read_qc/post-QC_reports
            else

                jobID=$(printf "%02d" $module)
                jobID=${jobID}_${ID}_long

                export moduleDir=../Duke_long/${moduleDir}
                logDir=${moduleDir}/${ID}
                mkdir -p $logDir
                mkdir -p ../Duke_long/${clean_readDir}
                mkdir -p ../Duke_long/${moduleDir}/pre-QC_reports
                mkdir -p ../Duke_long/${moduleDir}/post-QC_reports


                queued=$(job_lookup $jobID)
                # If this job is not already in the queue
                if [ -z "$queued"]; then
                    
                    # If the dependent job is complete
                    if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                        
                        count=$((count + 1))
                        H2 "${ID}"
                        H3 "Long-read QC"
                        LONG_QC_JOB=$(sbatch \
                            $long_read_qc_opts \
                            --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${ROOT}/pipelineScripts/scripts/long_scripts/${script})
                        echo $LONG_QC_JOB
                        echo "[ Log file ] -> ${logDir}/${ID}.${LONG_QC_JOB##* }.log"
                    
                    # If the dependent job is NOT complete
                    else

                        H3 "Long-read QC"
                        LONG_QC_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                            $long_read_qc_opts \
                            --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${ROOT}/pipelineScripts/scripts/long_scripts/${script})
                        echo $LONG_QC_JOB
                        echo "[ Log file ] -> ${logDir}/${ID}.${LONG_QC_JOB##* }.log"

                    fi # If the dependent job is complete
                
                # If this job is already in the queue
                else

                    LONG_QC_JOB=$queued

                fi # If this job is not already in the queue



            fi # If *_ont.fastq is already in Duke_long/clean_reads

        fi # If *_ont.fastq is already in HYBRID/clean_reads
        ((module++))
        if [[ "$1" = "read_qc" ]]; then continue; fi

    ##- Assembly
        
        script=assembly.sh
        moduleDir="${script%.sh}"
        export assemblyDir=$moduleDir

        assembly_OUT=${moduleDir}/${ID}/assembly.fasta

        # If assembly.fasta is already in assembly/*
        if [[ -s "$assembly_OUT" ]]; then

            ASSEMBLY_JOB=COMPLETE

        # If assembly.fasta is NOT already in assembly/*
        else

            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            # If this job is not already in the queue
            if [ -z "$queued" ]; then
                
                # If the read_QC module for both long and short has already been run
                if [ "$LONG_QC_JOB" = "COMPLETE" ] && [ "$SHORT_QC_JOB" = "COMPLETE" ]; then
                    
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Assembly"
                    # Submit a job with no dependencies
                    ASSEMBLY_JOB=$(sbatch \
                    $asm_opts \
                    --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $ASSEMBLY_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ASSEMBLY_JOB##* }.log"

                # If the long-read QC job is complete, but the short-read QC job is in progress    
                elif [ "$LONG_QC_JOB" = "COMPLETE" ] && [ "$SHORT_QC_JOB" != "COMPLETE" ]; then

                    H3 "Assembly"
                    # Submit a job dependent on the completion of the short-read QC
                    ASSEMBLY_JOB=$(sbatch --dependency=afterok:${SHORT_QC_JOB##* } \
                    $asm_opts \
                    --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $ASSEMBLY_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ASSEMBLY_JOB##* }.log"

                # If the long-read QC job is in progress, but the short-read QC job is complete  
                elif [ "$LONG_QC_JOB" != "COMPLETE" ] && [ "$SHORT_QC_JOB" = "COMPLETE" ]; then  

                    H3 "Assembly"
                    # Submit a job dependent on the completion of the long-read QC
                    ASSEMBLY_JOB=$(sbatch --dependency=afterok:${LONG_QC_JOB##* } \
                    $asm_opts \
                    --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $ASSEMBLY_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ASSEMBLY_JOB##* }.log"

                # If both the long-read and short-read QC jobs are in progress  
                else

                    H3 "Assembly"
                    # Submit a job dependent on the completion of both the long-read and short-read QC jobs
                    ASSEMBLY_JOB=$(sbatch --dependency=afterok:${LONG_QC_JOB##* }:${SHORT_QC_JOB##* } \
                    $asm_opts \
                    --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $ASSEMBLY_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ASSEMBLY_JOB##* }.log"

               fi # If the read_QC module for both long and short has already been run

            # If this job is already in the queue
            else

                ASSEMBLY_JOB=$queued

            fi # If this job is not already in the queue


        fi # If assembly.fasta is already in assembly/*
        
        ((module++))

        # if [[ $count -gt 30 ]]; then exit 0; fi
        if [[ "$1" = "assembly" ]]; then continue; fi

    ##- Evaluation

        script=evaluation.sh
        moduleDir="${script%.sh}"
        export evaluationDir=$moduleDir

        DEPENDENT_JOB=$ASSEMBLY_JOB

        eval_OUT=${moduleDir}/${ID}_final_assembly.fasta

        if [[ -s "$eval_OUT" ]]; then
            EVALUATION_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Evaluation"
                    EVALUATION_JOB=$(sbatch \
                        $eval_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $EVALUATION_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${EVALUATION_JOB##* }.log"
                else
                    H3 "Evaluation"
                    EVALUATION_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $eval_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $EVALUATION_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${EVALUATION_JOB##* }.log"
                fi
             else
                EVALUATION_JOB=$queued
          fi
       fi
       ((module++))

        if [[ "$1" = "evaluation" ]]; then continue; fi
    
    ##- Kraken2
        script=kraken2.sh
        moduleDir="${script%.sh}"

        DEPENDENT_JOB=$EVALUATION_JOB
        bracken_OUT=bracken/assembly/${ID}_assembly.bracken.out

        if [[ -s "$bracken_OUT" ]]; then
            KRAKEN_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/${ID}
            mkdir -p $logDir
            mkdir -p bracken

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Kraken2"
                    KRAKEN_JOB=$(sbatch \
                        $k2_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $KRAKEN_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${KRAKEN_JOB##* }.log"
                else
                    H3 "Kraken2"
                    KRAKEN_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $k2_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $KRAKEN_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${KRAKEN_JOB##* }.log"
                fi
            else
                KRAKEN_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "kraken2" ]]; then continue; fi

    ##- Binning
        
        script=binning.sh
        moduleDir="${script%.sh}"
        export binningDir=$moduleDir

        DEPENDENT_JOB=$EVALUATION_JOB

        metabat2_OUT=${moduleDir}/${ID}/metabat2_bins
        maxbin2_OUT=${moduleDir}/${ID}/maxbin2_bins
        concoct_OUT=${moduleDir}/${ID}/concoct_bins

        if [[ -d "$concoct_OUT" ]]; then
            BINNING_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then

                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Binning"
                    BINNING_JOB=$(sbatch \
                        $bin_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $BINNING_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${BINNING_JOB##* }.log"
                else
                    H3 "Binning"
                    BINNING_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $bin_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                    echo $BINNING_JOB
                    echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${BINNING_JOB##* }.log"
                fi
            else
                BINNING_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "binning" ]]; then continue; fi

    ##- Bin Refinement

        script=refine_bins.sh
        moduleDir="${script%.sh}"
        export refinedDir=$moduleDir

        export max_completion
        export min_contam

        
        refine_OUT=${moduleDir}/${ID}/metawrap_${max_completion}_${min_contam}_bins.stats

        if [[ -s "$refine_OUT" ]]; then
            REFINE_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$BINNING_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Bin Refinement"
                    REFINE_JOB=$(sbatch \
                        $refine_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $REFINE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${REFINE_JOB##* }.log"
                else
                    H3 "Bin Refinement"
                    REFINE_JOB=$(sbatch --dependency=afterok:${BINNING_JOB##* } \
                        $refine_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $REFINE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${REFINE_JOB##* }.log"
               fi
            else
                REFINE_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "refine" ]]; then continue; fi

    ##- Bin Reassembly

        script=reassemble_bins.sh
        moduleDir="${script%.sh}"
        export reassemDir=$moduleDir

        DEPENDENT_JOB=$REFINE_JOB
        reassembly_OUT=${moduleDir}/${ID}/reassembled_bins.stats

        if [[ -s "$reassembly_OUT" ]]; then
            REASSEMBLY_JOB=COMPLETE
        else
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Bin Reassembly"
                    REASSEMBLY_JOB=$(sbatch \
                        $reassem_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $REASSEMBLY_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${REASSEMBLY_JOB##* }.log"
                else
                    H3 "Bin Reassembly"
                    REASSEMBLY_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $reassem_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $REASSEMBLY_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${REASSEMBLY_JOB##* }.log"
               fi
            else
                REASSEMBLY_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "reassembly" ]]; then continue; fi

    ##- Bin Classification

        script=classify_bins.sh
        moduleDir="${script%.sh}"

        DEPENDENT_JOB=$REASSEMBLY_JOB

        classify_OUT=${moduleDir}/${ID}/bin_taxonomy.tab

        if [[ -s "$classify_OUT" ]]; then
            CLASSIFY_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Bin Classification"
                    CLASSIFY_JOB=$(sbatch \
                        $classify_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $CLASSIFY_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${CLASSIFY_JOB##* }.log"
                else
                    H3 "Bin Classification"
                    CLASSIFY_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $classify_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $CLASSIFY_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${CLASSIFY_JOB##* }.log"
               fi
            else
                CLASSIFY_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "classify" ]]; then continue; fi

    ##- Functional Annotation

        script=annotate_bins.sh
        moduleDir="${script%.sh}"

        DEPENDENT_JOB=$REASSEMBLY_JOB
        annotate_OUT=${moduleDir}/${ID}/COMPLETE

        if [[ -f "$annotate_OUT" ]]; then
            ANNOTATE_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "Functional Annotation"
                    ANNOTATE_JOB=$(sbatch \
                        $annotate_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $ANNOTATE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ANNOTATE_JOB##* }.log"
                else
                    H3 "Functional Annotation"
                    ANNOTATE_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
                        $annotate_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $ANNOTATE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ANNOTATE_JOB##* }.log"
               fi
            else
                ANNOTATE_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "annotation" ]]; then continue; fi

    ##- AMR Identification

        script=amr_detection.sh
        moduleDir="${script%.sh}"

        mkdir -p ${moduleDir}/AMR
        mkdir -p ${moduleDir}/RGI
        amr_OUT=${moduleDir}/RGI/${ID}.rgi.txt

        if [[ -s "$amr_OUT" ]]; then
            AMR_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$EVALUATION_JOB" = "COMPLETE" ]]; then
                    count=$((count + 1))
                    H2 "${ID}"
                    H3 "AMR Detection"
                    AMR_JOB=$(sbatch \
                        $amr_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $AMR_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${AMR_JOB##* }.log"
                else
                    H3 "AMR Detection"
                    AMR_JOB=$(sbatch --dependency=afterok:${EVALUATION_JOB##* } \
                        $amr_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $AMR_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${AMR_JOB##* }.log"
               fi
            else
                AMR_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "amr_search" ]]; then continue; fi

    ##- Assembly Gene Profiling
    
        inputType=asm
        script=gene_profiling.sh
        moduleDir="${inputType}_${script%.sh}"
        DEPENDENT_JOB=$EVALUATION_JOB

        mkdir -p ${moduleDir}/
        asm_OUT=${moduleDir}/$ID/${ID}_gene_annotations.tsv

        if [[ -s "$asm_OUT" ]]; then
            ASM_PROFILE_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export inputType
            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    if [[ "$AMR_JOB" = "COMPLETE" ]]; then count=$((count + 1)); H2 "${ID}"; fi
                    
                    H3 "Assembly Gene Profiling"
                    ASM_PROFILE_JOB=$(sbatch \
                        $asm_profile_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $ASM_PROFILE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ASM_PROFILE_JOB##* }.log"
                else
                    H3 "Assembly Gene Profiling"
                    ASM_PROFILE_JOB=$(sbatch --dependency=afterok:${EVALUATION_JOB##* } \
                        $asm_profile_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $ASM_PROFILE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${ASM_PROFILE_JOB##* }.log"
               fi
            else
                ASM_PROFILE_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "asm_profile" ]]; then continue; fi
    
    ##- Bin Gene Profiling
    
        inputType=bin
        script=gene_profiling.sh
        moduleDir="${inputType}_${script%.sh}"

        DEPENDENT_JOB=$EVALUATION_JOB

        mkdir -p ${moduleDir}/
        bin_OUT=${moduleDir}/$ID/${ID}_gene_annotations.tsv

        if [[ -s "$bin_OUT" ]]; then
            BIN_PROFILE_JOB=COMPLETE
        else
            
            jobID=$(printf "%02d" $module)
            jobID=${jobID}_${ID}_${readType}

            export inputType
            export moduleDir
            logDir=${moduleDir}/logs
            mkdir -p $logDir

            queued=$(job_lookup $jobID)
            if [ -z "$queued" ]; then
                if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
                    if [[ "$ASM_PROFILE_JOB" = "COMPLETE" ]]; then count=$((count + 1)); H2 "${ID}"; fi
                    H3 "Bin Gene Profiling"
                    BIN_PROFILE_JOB=$(sbatch \
                        $asm_profile_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $BIN_PROFILE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${BIN_PROFILE_JOB##* }.log"
                else
                    H3 "Bin Gene Profiling"
                    BIN_PROFILE_JOB=$(sbatch --dependency=afterok:${EVALUATION_JOB##* } \
                        $asm_profile_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
                   echo $BIN_PROFILE_JOB
                   echo "[ Log file ] -> HYBRID/${logDir}/${ID}.${BIN_PROFILE_JOB##* }.log"
               fi
            else
                BIN_PROFILE_JOB=$queued
            fi
        fi
        ((module++))

        if [[ "$1" = "bin_profile" ]]; then continue; fi

    # if [[ $count -ge 1 ]]; then exit 0; fi
done



