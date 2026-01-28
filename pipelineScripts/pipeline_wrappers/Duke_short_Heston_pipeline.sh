#!/bin/bash

    #-------------------------------------------------------------------------------------------------------------------
    # This script is used to begin the HCT ARG Assembly pipeline for the Illumina (short) read samples.

    # It is important to note that it has been designed for a specific working directory. Therefore, the reproduction of the results will require small modifications of the script or the adaptation of your working directory.

    # Created on Nov 6, 2024

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 3


    # Required tools:
        # 1. metaWRAP metagenomic wrapper suite (https://github.com/bxlab/metaWRAP)
        # 2. BBMap bioinformatic tools (https://github.com/BioInfoTools/BBMap)

    # This pipeline requires metaWRAP (https://github.com/bxlab/metaWRAP) for modules 00-01, 03-08. 

    # Ensure that the metaWRAP config file include paths to the necessary databases needed for this pipeline.
    # Instructions found here: https://github.com/bxlab/metaWRAP/blob/master/installation/database_installation.md

    # Database      Size     Used in module
    # -----------------------------------------
    # CheckM        1.4Gb    binning, bin_refinement, reassemble_bins
    # Kraken2       125Gb    kraken2
    # NCBI_nt       71Gb     classify_bins
    # NCBI_tax      283Mb    classify_bins
    # Indexed hg38  20Gb     read_qc

    # Module 02 requires bbmap

    # Module 09 requires the tools AMRFinder+ (https://github.com/ncbi/amr) and RGI (https://github.com/arpcard/rgi) with the required databases.

    # The sampleList is a text file of the raw sequence sample names in the following format:
    # #SampleID
    # D20248PRE_GGACTCCT-ACTGCATA_S221_L004
    # D20248D1_CTCTCTAC-ACTGCATA_S134_L003

    # This pipeline was originally run on Red Hat Enterprise Linux 9.2 (Plow) using the Slurm Workload Manager.

    #-------------------------------------------------------------------------------------------------------------------

cohort=Duke
read=short
dataset=${cohort}_${read}

source pipelineScripts/configs/${dataset}-read.config
export bashrc
export pipelineConfig=${ROOT}/pipelineScripts/configs/${dataset}-read.config
export config_file=$(which config-metawrap)

if [[ ! -d $ROOT/${dataset} ]]; then mkdir -p $ROOT/${dataset}; fi
cd $ROOT/${dataset}

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

        if [ ! -z "$2" ]; then moduleDir=$moduleDir/$2; fi

        if [[ ! -d ${moduleDir} ]]; then mkdir -p $moduleDir; fi
        export moduleDir
        
        logDir=${moduleDir}/logs
        if [[ ! -d ${logDir} ]]; then mkdir -p $logDir; fi

        Complete_tag=(${moduleDir}/COMPLETE/${ID})

        jobID="${script%%_*}"
        jobID=${jobID}_${ID}_${dataset}
        if [ ! -z "$2" ]; then jobID=${jobID}_$2; fi

    }
    first_ID() {
        if [[ $module == 0 ]]; then 
            H2 "${ID}"
            ((module++))
            ((count++))
        fi    
    }
    run_module() {
        
        # If the completion file exists
        all_exist=true

        for file in "${Complete_tag[@]}"; do
          if [[ ! -e "$file" ]]; then
            echo $file
            all_exist=false
            break
          fi
        done

        if $all_exist; then
            
            Current_Job=COMPLETE

        else

            queued=$(job_lookup $jobID)
            # If the job is not queued
            if [ -z "$queued" ]; then

                ready_to_run=true
                dependencies=""
                for job in "${DEPENDENT_JOB[@]}"; do
                    if [[ $job != "COMPLETE" ]]; then
                        ready_to_run=false
                        if [[ -z "$dependencies" ]]; then
                            dependencies="${job}"
                        else
                            dependencies+=",${job}"
                        fi
                    fi
                done
                echo $dependencies

                # If the dependent job is complete
                if $ready_to_run; then

                    first_ID
                    H3 "$header3"
                    Current_Job=$(sbatch \
                        $hpc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${moduleScriptPath}/${script})
                    echo $Current_Job
                    echo "[ Log file ] -> ${dataset}/${logDir}/${ID}.${Current_Job##* }.log"

                else
                    
                    H3 "$header3"
                    # echo -e "sbatch --dependency=afterok:${dependencies} \
                    #     $hpc_opts \
                    #     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${moduleScriptPath}/${script}"
                    Current_Job=$(sbatch --dependency=afterok:${dependencies} \
                        $hpc_opts \
                        --job-name=${jobID} -o ${logDir}/${ID}.%A.log ${moduleScriptPath}/${script})
                    echo $Current_Job
                    echo "[ Log file ] -> ${dataset}/${logDir}/${ID}.${Current_Job##* }.log"

               fi 

            else
                Current_Job=$queued
           fi # if [ -z "$queued" ]

        fi     
    }

# Set up
    if [[ ! -f "LOGs/${dataset}_pipeline_$version.out" ]]; then
        H3 "Usage"
        echo "export version=$version"
        echo "nohup sh ./pipelineScripts/${dataset}_pipeline.sh >> ${dataset}/LOGs/${dataset}_pipeline_$version.out 2>&1 &"
        echo 
        comment "[ Raw sequence directory ]: ${seqPath}"
    fi
        

    mkdir -p LOGs
    export seqPath
    export readType
    export totalSamples=`tail -n +2 $sampleList | wc -l`

    # Initialize count variable
    count=0

# Run pipeline
    for s in $(tail -n +2 $sampleList); do
        export s
        export ID="${s%%_*}"
        export raw_readDir=0.0_raw_reads
        export pre_qcDir=0.1_pre_qc
        export dedup_Dir=0.2_deduplication
        export trimmed_Dir=0.3_sequence_trim
        export clean_readDir=0.4_host_decontamination
        export krakenDir=2.1_kraken2
        export brackenDir=2.2_bracken
        export shortbredDir=5.3_shortbred

        module=0
        if [[ $ID == "D21309D98" || $ID == "BMT174PRE" ]]; then continue; fi
        # if [[ $count -ge 5 ]]; then continue; fi

        ##- 0.1 Pre-QC

            module_setup 0.1_pre_qc.sh

            # Module inputs -------------
            header3="Pre-QC"
            DEPENDENT_JOB=(COMPLETE)
            hpc_opts=$pre_qc_opts
            pipeline_tag=pre_qc
            # ---------------------------


            # Module specific -----------
            if [[ ! -d ${raw_readDir} ]]; then mkdir -p $raw_readDir; fi
            export R1_ext
            export R2_ext
            # Complete_tag=()
            # Complete_tag+=(${raw_readDir}/${ID}_1.fastq)
            # Complete_tag+=(${raw_readDir}/${ID}_2.fastq)
            # Complete_tag+=(${pre_qcDir}/${ID}_1_fastqc.html)
            # Complete_tag+=(${pre_qcDir}/${ID}_1_fastqc.zip)
            # Complete_tag+=(${pre_qcDir}/${ID}_2_fastqc.html)
            # Complete_tag+=(${pre_qcDir}/${ID}_2_fastqc.zip)
            # Complete_tag+=(${clean_readDir}/${ID}_1.fastq.gz)
            # Complete_tag+=(${clean_readDir}/${ID}_2.fastq.gz)
            # ---------------------------
            
            run_module
            PRE_QC_JOB=$Current_Job        
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


        ##- 0.2 Deduplication

            module_setup 0.2_deduplication.sh

            # Module inputs -------------
            header3="Deduplication"
            DEPENDENT_JOB=(${PRE_QC_JOB##* })
            hpc_opts=$dedup_opts
            pipeline_tag=dedup
            # ---------------------------


            # Module specific -----------
            # Complete_tag=()
            # Complete_tag+=(${moduleDir}/${ID}_stats.log)
            # Complete_tag+=(${moduleDir}/QC_report/${ID}_deduped_R1_fastqc.html)
            # Complete_tag+=(${moduleDir}/QC_report/${ID}_deduped_R1_fastqc.zip)
            # Complete_tag+=(${moduleDir}/QC_report/${ID}_deduped_R2_fastqc.html)
            # Complete_tag+=(${moduleDir}/QC_report/${ID}_deduped_R2_fastqc.zip)
            # Complete_tag+=(${clean_readDir}/${ID}_1.fastq.gz)
            # Complete_tag+=(${clean_readDir}/${ID}_2.fastq.gz)
            # ---------------------------
            
            run_module
            DEDUP_JOB=$Current_Job        
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


        ##- 0.3 Sequence Trimming

            module_setup 0.3_sequence_trim.sh

            # Module inputs -------------
            header3="Sequence Trimming"
            DEPENDENT_JOB=(${DEDUP_JOB##* })
            hpc_opts=$trim_opts
            pipeline_tag=trim
            # ---------------------------


            # Module specific -----------
            # ---------------------------
            
            run_module
            TRIM_JOB=$Current_Job        
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


        ##- 0.4 Host Decontamination

            module_setup 0.4_host_decontamination.sh

            # Module inputs -------------
            header3="Host Decontamination"
            DEPENDENT_JOB=(${TRIM_JOB##* })
            hpc_opts=$decontam_opts
            pipeline_tag=decontam
            # ---------------------------


            # Module specific -----------
            # ---------------------------
            
            run_module
            DECONTAM_JOB=$Current_Job    
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


        ##- 2.1 Kraken2

            module_setup 2.1_kraken2.sh

            # Module specific -------------------
            header3="Kraken2"
            DEPENDENT_JOB=(${DECONTAM_JOB##* })
            hpc_opts=$k2_opts
            pipeline_tag=kraken2
            # -----------------------------------

            # Module specific actions -----------
            if [[ ! -d $brackenDir ]]; then mkdir -p $brackenDir; fi
            # -----------------------------------

            
            run_module 
            KRAKEN_JOB=$Current_Job        
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


        ##- 5.3 ShortBRED AMR Identification
            
            module_setup 5.3_shortbred.sh
            
            # Module specific -------------------
            header3="ShortBRED AMR Identification"
            DEPENDENT_JOB=(${DECONTAM_JOB##* })
            hpc_opts=$shortbred_opts
            pipeline_tag=shortbred
            # -----------------------------------
            
            # Module specific actions -----------
            # -----------------------------------
            
            
            run_module 
            ShortBRED_JOB=$Current_Job       
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


        ##- 6.1 Discard intermediate files
            
            module_setup COMPLETE.sh
            
            # Module specific -------------------
            header3="Discard intermediate files"
            DEPENDENT_JOB=(${PRE_QC_JOB##* } ${DEDUP_JOB##* } ${TRIM_JOB##* } ${DECONTAM_JOB##* } ${KRAKEN_JOB##* } ${ShortBRED_JOB##* })
            hpc_opts=$discard_opts
            pipeline_tag=discard
            # -----------------------------------
            
            # Module specific actions -----------
            Complete_tag=(COMPLETE/${ID})
            jobID=${ID}_cleanup
            # -----------------------------------
            
            
            run_module 
            DISCARD_JOB=$Current_Job       
            if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


    done
    H1 "Pipeline Complete!"


