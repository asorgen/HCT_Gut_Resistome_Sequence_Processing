#!/bin/bash

# Description
    # This script is used to begin the HCT ARG Assembly pipeline for the Oxford Nanopore (long) read samples.

    # It is important to note that it has been designed for a specific working directory. Therefore, the reproduction
    # of the results will require small modifications of the script or the adaptation of your working directory.

    # Created on Nov 6, 2024; Updated March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 4

    # Required tools:
        # 1. metaWRAP metagenomic wrapper suite (https://github.com/bxlab/metaWRAP)
        # 2. Porechop ONT adapter trimmer (https://github.com/rrwick/Porechop)
        # 3. NanoFilt read quality/length filter (https://github.com/wdecoster/nanofilt)
        # 4. minimap2 sequence aligner (https://github.com/lh3/minimap2)
        # 5. samtools sequence manipulation (https://github.com/samtools/samtools)
        # 6. Flye de novo PacBio/ONT genome assembler (https://github.com/mikolmogorov/Flye)
        # 7. AMRFinder+ resistance gene identifier (https://github.com/ncbi/amr)
        # 8. RGI resistance gene identifier (https://github.com/arpcard/rgi)
        # 9. Prodigal ORF prediction (https://github.com/hyattpd/Prodigal)
        # 10. GTDBtk genome taxonomy (https://github.com/Ecogenomics/GTDBTk)
        # 11. Bakta functional annotation (https://github.com/oschwengers/bakta)
        # 12. Kraken2 + Bracken taxonomic classification

    # ONT-specific module notes:
        # 0.2 Deduplication: Not included (ligation library prep does not introduce PCR duplicates)
        # 0.3 Sequence trim: Porechop (adapter removal) + NanoFilt (quality/length filter)
        # 0.4 Host decontam: minimap2 -ax map-ont to GRCh38 (replaces Bowtie2)
        # 1.1 Assembly: metaFlye --nano-hq (replaces metaSPAdes)
        # 2.1 Kraken2: single-end ONT reads (no --paired; assembly classification unchanged)
        # 3.1 Binning: minimap2 -x map-ont alignment for coverage (replaces BWA-MEM)
        # 3.3 Reassembly: Flye (replaces SPAdes)
        # 5.2 NT AMR (bins): Not included in long-read pipeline
        # 5.4 RGI BWT: single-end ONT FASTQ (no --read_two), KMA aligner
        # 5.6 AA AMR (bins): Not included in long-read pipeline

    # The sampleList is a text file of the sample names in the following format:
    # #SampleID
    # D20248PRE
    # D20248D1

    # This pipeline was originally run on Red Hat Enterprise Linux 9.2 (Plow) using the Slurm Workload Manager.


# Set pipeline 
    cohort=Duke
    read=long
    dataset=${cohort}_${read}

# Source configs and functions
    source pipelineScripts/configs/${dataset}-read.config
    source pipelineScripts/configs/functions.sh
    export bashrc
    export pipelineConfig=${ROOT}/pipelineScripts/configs/${dataset}-read.config
    export config_file=$(which config-metawrap)
    export module_functions
    export print_functions


# Set up
    if [[ ! -d $datasetDir ]]; then mkdir -p $datasetDir; fi
    cd $datasetDir

    if [[ ! -f "LOGs/${dataset}_pipeline_$version.out" ]]; then
        H3 "Usage"
        echo "export version=$version"
        echo "nohup sh ./pipelineScripts/${dataset}_pipeline.sh >> ${dataset}/LOGs/${dataset}_pipeline_$version.out 2>&1 &"
        echo 
        comment "[ Raw sequence directory ]: ${seqPath}"
    fi
        
    mkdir -p LOGs

    # Capture stop-tag argument so it remains accessible inside run_module_step
    pipeline_stop_tag="$1"

    # Export all module directories so sbatch jobs inherit them
    # (replaces per-module `export varName` lines in each block)
    export raw_readDir pre_qcDir trimmed_Dir clean_readDir
    export assemblyDir evaluationDir krakenDir brackenDir
    export binningDir refinedbinDir reassemDir rgi_bwt_dir
    export GRCh38_FASTA

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
        if [[ $count -ge 1 ]]; then continue; fi

        ##- 0.0 Copy raw Nanopore file
            if $run_cp_raw; then
                hpc_opts="--partition=Orion --nodes=1 --cpus-per-task=4 --mem=8GB --time=24:00:00 --mail-user=${email} --mail-type=FAIL"
                module_step_run "0.0_raw_reads.sh" "Copy raw Nanopore files" \
                    (COMPLETE) "$hpc_opts" "cp_raw" "RAW_READS_JOB" || continue

                # module_setup 0.0_raw_reads.sh

                # # Module inputs -------------
                # header3="Copy raw Nanopore files"
                # DEPENDENT_JOB=(COMPLETE)
                # hpc_opts="--partition=Orion --nodes=1 --cpus-per-task=4 --mem=8GB --time=24:00:00 --mail-user=${email} --mail-type=FAIL"
                # pipeline_tag=cp_raw
                # if [[ ! -d ${raw_readDir} ]]; then mkdir -p $raw_readDir; fi
                # export raw_readDir
                # # Skip copy if pre-QC is already done
                # Complete_tag=(${raw_readDir}/${ID}.fastq.gz)
                # # -----------------------------------

                # run_module
                # RAW_READS_JOB=$Current_Job
                # if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
            fi

        ##- 0.1 Pre-QC
            if $run_pre_qc; then
                run_module_step "0.1_pre_qc.sh" "Pre-QC" \
                    "${RAW_READS_JOB##* }" "$pre_qc_opts" "pre_qc" "PRE_QC_JOB" || continue
            fi

        ##- 0.3 Sequence trim (Porechop adapter trimming + NanoFilt quality/length filtering)
            if $run_trim; then
                run_module_step "0.3_sequence_trim.sh" "Sequence trim (Porechop + NanoFilt)" \
                    "${PRE_QC_JOB##* }" "$trim_opts" "trim" "TRIM_JOB" || continue
            fi

        ##- 0.4  Host decontamination (minimap2 map-ont + samtools)
            if $run_decontam; then
                run_module_step "0.4_host_decontamination.sh" "Host decontamination (minimap2 map-ont)" \
                    "${TRIM_JOB##* }" "$decontam_opts" "decontam" "DECONTAM_JOB" || continue
            fi


        ##- 1.1 Assembly (metaFlye)
            if $run_asm; then
                run_module_step "1.1_assembly.sh" "Assembly (metaFlye)" \
                    "${DECONTAM_JOB##* }" "$asm_opts" "assembly" "ASSEMBLY_JOB" || continue
            fi

        ##- 1.2 Evaluation (filter contigs <${min_contig_len}bp)
            if $run_asm; then
                run_module_step "1.2_evaluation.sh" "Evaluation (filter short contigs)" \
                    "${ASSEMBLY_JOB##* }" "$eval_opts" "eval" "EVALUATION_JOB" || continue
            fi

    done
    H1 "Pipeline Complete!"


