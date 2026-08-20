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


# Set pipeline
    # The cohort is no longer baked in here. The config declares cohort/read/dataset
    # along with every path, so adding a cohort means writing a config, not copying
    # this file. Accepts either a bare config name or a path:
    #   ./run_short_pipeline.sh Heston_short
    #   ./run_short_pipeline.sh Heston_short trim
    #   ./run_short_pipeline.sh pipelineScripts/configs/Heston_short-read.config

    CONFIG_ARG=$1
    if [[ -z "$CONFIG_ARG" ]]; then
        echo "Usage: run_short_pipeline.sh <config> [stop_at]"
        echo "  <config>   config name (e.g. Heston_short) or path to a *-read.config"
        echo "  [stop_at]  optional pipeline tag to halt after (e.g. trim, kraken, rgi)"
        echo
        echo "Available short-read configs:"
        ls pipelineScripts/configs/*_short-read.config 2>/dev/null \
            | xargs -n1 basename 2>/dev/null | sed 's/-read\.config$//; s/^/  /'
        exit 1
    fi

    # The stop_at tag is compared against "$1" throughout the body, so shift the
    # config off and leave stop_at sitting in $1.
    shift

    if [[ -f "$CONFIG_ARG" ]]; then
        CONFIG_FILE=$CONFIG_ARG
    else
        CONFIG_FILE=pipelineScripts/configs/${CONFIG_ARG}-read.config
    fi

    if [[ ! -f "$CONFIG_FILE" ]]; then
        echo "No such config: ${CONFIG_FILE}"
        echo "Run without arguments to list the available configs."
        exit 1
    fi

# Source configs and functions
    source "$CONFIG_FILE"
    source pipelineScripts/configs/functions.sh
    export bashrc
    export pipelineConfig=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")
    export config_file=$(which config-metawrap)
    export module_functions
    export print_functions

# Fail-safe module flags
    # An undefined flag makes `if $run_x` expand to an empty command, which bash
    # treats as TRUE — so a config that simply omits a module would silently RUN it.
    # Defaulting every flag to false after sourcing makes omission mean "skip".
    for _flag in pre_qc dedup trim decontam asm eval k2 m4 bin refine reassem \
                 classify annotate amr_nt_asm amr_nt_bin shortbred rgi_bwt \
                 amr_aa_asm amr_aa_bin waafle humann clean_up; do
        eval ": \${run_${_flag}:=false}"
    done
    unset _flag

# Optional sample-list override
    # Lets a run cover a subset of the cohort without editing the config or
    # duplicating it. Everything else (datasetDir, COMPLETE flags, module dirs)
    # still comes from the config, so batches accumulate into the same run.
    # Needed because Orion's QOS caps submissions at 2048 and the wrapper submits
    # one job per module per sample, so a large cohort must go in batches.
    if [[ -n "${SAMPLE_LIST_OVERRIDE:-}" ]]; then
        if [[ ! -f "$SAMPLE_LIST_OVERRIDE" ]]; then
            echo "SAMPLE_LIST_OVERRIDE set but not a file: ${SAMPLE_LIST_OVERRIDE}"
            exit 1
        fi
        sampleList=$SAMPLE_LIST_OVERRIDE
        echo "Sample list overridden: ${sampleList}"
    fi

# Optional per-invocation submission cap
    # 0 = unlimited, which is the historical behaviour, so omitting this changes
    # nothing for Duke/UNC. `count` is incremented by first_ID(), which run_module
    # calls only when it actually submits a sample's first module — samples whose
    # modules are all COMPLETE short-circuit before that. So the cap is
    # self-advancing: re-running submits the NEXT N samples still needing work.
    #
    # NOTE: this caps submissions PER INVOCATION, not total in-flight. Re-running
    # while the previous N are still queued will submit N more, because run_module
    # sees them as queued rather than needing work. Wait for the queue to drain
    # between invocations (see run_short_pipeline_throttled.sh).
    : "${MAX_SAMPLES_PER_RUN:=0}"
    if [[ ! $MAX_SAMPLES_PER_RUN =~ ^[0-9]+$ ]]; then
        echo "MAX_SAMPLES_PER_RUN must be a non-negative integer: ${MAX_SAMPLES_PER_RUN}"
        exit 1
    fi
    [[ $MAX_SAMPLES_PER_RUN -gt 0 ]] && \
        echo "Submission cap: ${MAX_SAMPLES_PER_RUN} samples this invocation"

    if [[ -z "$cohort" || -z "$read" || -z "$dataset" ]]; then
        echo "Config ${CONFIG_FILE} does not set cohort/read/dataset."
        exit 1
    fi

    if [[ "$read" != "short" ]]; then
        echo "Config ${CONFIG_FILE} declares read=${read}; this wrapper is short-read only."
        exit 1
    fi


# Set up
    if [[ ! -d $datasetDir ]]; then mkdir -p $datasetDir; fi
    cd $datasetDir

    if [[ ! -f "LOGs/${dataset}_pipeline_$version.out" ]]; then
        H3 "Usage"
        echo "export version=$version"
        echo "nohup sh ./pipelineScripts/pipeline_wrappers/run_short_pipeline.sh ${dataset} >> ${dataset}/LOGs/${dataset}_pipeline_$version.out 2>&1 &"
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
        if echo "${exclude_ids}" | grep -qw "$ID"; then continue; fi
        if [[ $MAX_SAMPLES_PER_RUN -gt 0 && $count -ge $MAX_SAMPLES_PER_RUN ]]; then continue; fi

        ##- 0.1 Pre-QC
            if $run_pre_qc; then
                module_setup 0.1_pre_qc.sh

                # Module inputs -------------
                header3="Pre-QC"
                DEPENDENT_JOB=(COMPLETE)
                hpc_opts=$pre_qc_opts
                pipeline_tag=pre_qc
                if [[ ! -d ${raw_readDir} ]]; then mkdir -p $raw_readDir; fi
                export raw_readDir
                export pre_qcDir
                export R1_ext
                export R2_ext
                # ---------------------------
                
                run_module
                PRE_QC_JOB=$Current_Job        
                CLEAN_UP_DEP+=(${Current_Job##* })
                # echo ${CLEAN_UP_DEP[@]}
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
                export dedup_Dir
                # ---------------------------
                
                run_module
                DEDUP_JOB=$Current_Job        
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export trimmed_Dir
                # ---------------------------
                
                run_module
                TRIM_JOB=$Current_Job        
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export clean_readDir
                # ---------------------------
                
                run_module
                DECONTAM_JOB=$Current_Job    
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export assemblyDir
                # ---------------------------
                
                run_module
                ASSEMBLY_JOB=$Current_Job    
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export evaluationDir
                # ---------------------------
                
                run_module
                EVALUATION_JOB=$Current_Job    
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export krakenDir
                export brackenDir
                if [[ ! -d $brackenDir ]]; then mkdir -p $brackenDir; fi
                # -----------------------------------
                
                run_module 
                KRAKEN_JOB=$Current_Job        
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export metaphlanDir
                # ----------------------------
                
                run_module
                METAPHLAN_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export binningDir
                # ---------------------------
                
                run_module
                BINNING_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export refinedbinDir
                export min_completion; export max_contam
                # ---------------------------
                
                run_module
                REFINE_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export reassemDir
                # ---------------------------
                
                run_module
                REASSEMBLE_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export shortbredDir
                # -----------------------------------
                
                run_module 
                ShortBRED_JOB=$Current_Job       
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                export rgi_bwt_dir
                # -----------------------------------
                
                run_module 
                RGI_BWT_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                CLEAN_UP_DEP+=(${Current_Job##* })
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
                CLEAN_UP_DEP+=(${Current_Job##* })
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi
            

        ##- 5.7 WAAFLE — lateral gene transfer detection from assembled contigs
            if $run_waafle; then
                module_setup 5.7_waafle.sh

                # Module inputs -------------
                header3="WAAFLE LGT detection"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$waafle_opts
                pipeline_tag=waafle
                export evaluationDir
                export CHOCOPHLAN_BLAST_DB
                export CHOCOPHLAN_TAXONOMY
                export WAAFLE_ENV
                export waafleDir
                # ---------------------------

                run_module
                WAAFLE_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 5.8 HUMAnN — whole-community gene-family/pathway profiling
            if $run_humann; then
                module_setup 5.8_humann.sh

                # Module inputs -------------
                header3="HUMAnN gene-family/pathway profiling"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$humann_opts
                pipeline_tag=humann
                export clean_readDir
                export HUMANN_CHOCOPHLAN
                export HUMANN_UNIREF
                export humannDir
                # ---------------------------

                run_module
                HUMANN_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
                if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi


        ##- 6.1 Discard intermediate files
            if $run_clean_up; then
                module_setup COMPLETE.sh

                # Module specific -------------------
                header3="Discard intermediate files"
                DEPENDENT_JOB=(${CLEAN_UP_DEP[@]})
                # echo ${DEPENDENT_JOB[@]}
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


