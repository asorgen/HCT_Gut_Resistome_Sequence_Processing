#!/bin/bash

# Description
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


cohort=Duke
read=short
dataset=${cohort}_${read}

source pipelineScripts/configs/${dataset}-read.config
export bashrc
export pipelineConfig=${ROOT}/pipelineScripts/configs/${dataset}-read.config
export config_file=$(which config-metawrap)
export module_functions

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
        CLEAN_UP_DEP=()
        
        module=0
        # if [[ $ID == "D21309D98" ]]; then continue; fi
        if [[ $ID == "D21309D98" || $ID == "D13004PRE" ]]; then continue; fi
        if [[ $count -ge 50 ]]; then continue; fi

        ##- 0.1 Pre-QC
            if $run_pre_qc; then
                module_setup 0.1_pre_qc.sh

                # Module inputs -------------
                header3="Pre-QC"
                DEPENDENT_JOB=(COMPLETE)
                hpc_opts=$pre_qc_opts
                pipeline_tag=pre_qc
                export raw_readDir=0.0_raw_reads
                if [[ ! -d ${raw_readDir} ]]; then mkdir -p $raw_readDir; fi
                export pre_qcDir=$moduleDir
                export R1_ext
                export R2_ext
                # ---------------------------
                
                run_module
                PRE_QC_JOB=$Current_Job        
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 0.2 Deduplication
            if $run_dedup; then 
                module_setup 0.2_deduplication.sh

                # Module inputs -------------
                header3="Deduplication"
                DEPENDENT_JOB=(${PRE_QC_JOB##* })
                hpc_opts=$dedup_opts
                pipeline_tag=dedup
                export dedup_Dir=$moduleDir
                # ---------------------------
                
                run_module
                DEDUP_JOB=$Current_Job        
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 0.3 Sequence Trimming
            if $run_trim; then
                module_setup 0.3_sequence_trim.sh

                # Module inputs -------------
                header3="Sequence Trimming"
                DEPENDENT_JOB=(${DEDUP_JOB##* })
                hpc_opts=$trim_opts
                pipeline_tag=trim
                export trimmed_Dir=$moduleDir
                # ---------------------------
                
                run_module
                TRIM_JOB=$Current_Job        
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 0.4 Host Decontamination
            if $run_decontam; then
                module_setup 0.4_host_decontamination.sh

                # Module inputs -------------
                header3="Host Decontamination"
                DEPENDENT_JOB=(${TRIM_JOB##* })
                hpc_opts=$decontam_opts
                pipeline_tag=decontam
                export clean_readDir=$moduleDir
                # ---------------------------
                
                run_module
                DECONTAM_JOB=$Current_Job    
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 1.1 Assembly
            if $run_asm; then
                module_setup 1.1_assembly.sh

                # Module inputs -------------
                header3="Assembly"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$asm_opts
                pipeline_tag=assembly
                export assemblyDir=$moduleDir
                # ---------------------------
                
                run_module
                ASSEMBLY_JOB=$Current_Job    
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 1.2 Evaluation
            if $run_eval; then
                module_setup 1.2_evaluation.sh

                # Module inputs -------------
                header3="Evaluation"
                DEPENDENT_JOB=(${ASSEMBLY_JOB##* })
                hpc_opts=$eval_opts
                pipeline_tag=evaluation
                export evaluationDir=$moduleDir
                # ---------------------------
                
                run_module
                EVALUATION_JOB=$Current_Job    
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 2.1 Kraken2
            if $run_k2; then
                module_setup 2.1_kraken2.sh

                # Module specific -------------------
                header3="Kraken2"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$k2_opts
                pipeline_tag=kraken2
                export krakenDir=$moduleDir
                export brackenDir=2.2_bracken
                if [[ ! -d $brackenDir ]]; then mkdir -p $brackenDir; fi
                # -----------------------------------
                
                run_module 
                KRAKEN_JOB=$Current_Job        
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 2.3 MetaPhlAn4
            if $run_m4; then
                module_setup 2.3_metaphlan4.sh
                
                # Module specific -------------
                header3="MetaPhlAn4"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$m4_opts
                pipeline_tag=metaphlan4
                export metaphlanDir=$moduleDir
                # ----------------------------
                
                run_module
                METAPHLAN_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 3.1 Binning
            if $run_bin; then
                module_setup 3.1_binning.sh
                
                # Module specific -------------------
                header3="Binning"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$bin_opts
                pipeline_tag=binning
                export binningDir=$moduleDir
                # ---------------------------
                
                run_module
                BINNING_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 3.2 Refine Bins
            if $run_refine; then
                module_setup 3.2_refine_bins.sh
                
                # Module specific -------------------
                header3="Refine Bins"
                DEPENDENT_JOB=(${BINNING_JOB##* })
                hpc_opts=$refine_opts
                pipeline_tag=refine
                export refinedbinDir=$moduleDir
                export min_completion; export max_contam
                # ---------------------------
                
                run_module
                REFINE_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 3.3 Reassemble bins
            if $run_reassem; then
                module_setup 3.3_reassemble_bins.sh
                
                # Module specific -------------------
                header3="Reassemble bins"
                DEPENDENT_JOB=(${REFINE_JOB##* })
                hpc_opts=$reassem_opts
                pipeline_tag=reassemble
                export reassemDir=$moduleDir
                # ---------------------------
                
                run_module
                REASSEMBLE_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 4.1 Classify bins
            if $run_classify; then
                module_setup 4.1_classify_bins.sh
                
                # Module specific -------------------
                header3="Classify bins"
                DEPENDENT_JOB=(${REASSEMBLE_JOB##* })
                hpc_opts=$classify_opts
                pipeline_tag=classify
                #------------------------------------
                
                run_module
                ANNOTATE_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 4.2 Annotate bins
            if $run_annotate; then
                module_setup 4.2_annotate_bins.sh
                
                # Module specific -------------------
                header3="Annotate bins"
                DEPENDENT_JOB=(${REASSEMBLE_JOB##* })
                hpc_opts=$annotate_opts
                pipeline_tag=annotate
                # -----------------------------------
                
                run_module
                ANNOTATE_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.1 AMR detection from assembly nucleotide sequences
            if $run_amr_nt_asm; then
                module_setup 5.1_NT_amr_assembly.sh
                
                # Module specific -------------------
                header3="AMR detection from assembly nucleotide sequences"
                DEPENDENT_JOB=$EVALUATION_JOB
                hpc_opts=$amr_opts
                pipeline_tag=nt_asm_amr
                # -----------------------------------
                
                run_module 
                AMR_NT_ASM_JOB=$Current_Job       
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.2 AMR detection from MAG nucleotide sequences
            if $run_amr_nt_bin; then
                module_setup 5.2_NT_amr_bins.sh
                
                # Module specific -------------------
                header3="AMR detection from MAG nucleotide sequences"
                DEPENDENT_JOB=$REASSEMBLE_JOB
                hpc_opts=$amr_opts
                pipeline_tag=nt_bin_amr
                # -----------------------------------
                
                run_module 
                AMR_NT_BIN_JOB=$Current_Job       
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.3 ShortBRED AMR Identification
            if $run_shortbred; then
                module_setup 5.3_shortbred.sh
                
                # Module specific -------------------
                header3="ShortBRED AMR Identification"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$shortbred_opts
                pipeline_tag=shortbred
                export shortbredDir=$moduleDir
                # -----------------------------------
                
                run_module 
                ShortBRED_JOB=$Current_Job       
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.4 RGI BWT AMR Identification
            if $run_rgi_bwt; then
                module_setup 5.4_rgi_bwt.sh
                
                # Module specific -------------------
                header3="RGI BWT AMR Identification"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$rgi_bwt_opts
                pipeline_tag=rgi_bwt
                export aligner=kma
                export rgi_bwt_dir=$moduleDir
                # -----------------------------------
                
                run_module 
                RGI_BWT_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.5 AMR detection from assembled predicted genes
            if $run_amr_aa_asm; then
                module_setup 5.5_AA_amr_assembly.sh
                
                # Module specific -------------------
                header3="AMR detection from assembled predicted genes"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$asm_profile_opts
                pipeline_tag=aa_asm_amr
                # -----------------------------------
                
                run_module
                AMR_AA_ASM_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.6 AMR detection from assembled predicted genes
            if $run_amr_aa_bin; then
                module_setup 5.6_AA_amr_bins.sh
                
                # Module specific -------------------
                header3="AMR detection from assembled predicted genes"
                DEPENDENT_JOB=(${REASSEMBLE_JOB##* })
                hpc_opts=$asm_profile_opts
                pipeline_tag=aa_bin_amr
                # -----------------------------------
                
                run_module
                AMR_AA_BIN_JOB=$Current_Job
                CLEAN_UP_DEP+=($Current_Job)
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi
            

        ##- 6.1 Discard intermediate files
            if $run_clean_up; then
                module_setup COMPLETE.sh
                
                # Module specific -------------------
                header3="Discard intermediate files"
                DEPENDENT_JOB=($CLEAN_UP_DEP)
                hpc_opts=$discard_opts
                pipeline_tag=discard
                Complete_tag=(COMPLETE/${ID})
                jobID=${ID}_cleanup
                # -----------------------------------
                
                run_module 
                DISCARD_JOB=$Current_Job       
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


    done
    H1 "Pipeline Complete!"


