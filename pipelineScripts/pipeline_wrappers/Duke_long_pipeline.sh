#!/bin/bash

# Description
    #-------------------------------------------------------------------------------------------------------------------
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
        # 0.2 Deduplication: Not used (ligation library prep does not introduce PCR duplicates)
        # 0.3 Sequence trim: Porechop (adapter removal) + NanoFilt (quality/length filter)
        # 0.4 Host decontam: minimap2 -ax map-ont to GRCh38 (replaces Bowtie2)
        # 1.1 Assembly: metaFlye (replaces metaSPAdes)
        # 3.1 Binning: minimap2 -x map-ont alignment (replaces BWA-MEM)
        # 3.3 Reassembly: Flye (replaces SPAdes)
        # 5.4 RGI BWT: single-end ONT FASTQ (no --read_two)

    # The sampleList is a text file of the sample names in the following format:
    # #SampleID
    # D20248PRE
    # D20248D1

    # This pipeline was originally run on Red Hat Enterprise Linux 9.2 (Plow) using the Slurm Workload Manager.

    #-------------------------------------------------------------------------------------------------------------------

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
    export seqPath
    export readType
    export totalSamples=`tail -n +2 $sampleList | wc -l`

    # Initialize count variable
    count=0

# Run pipeline
    for ID in $(tail -n +2 $sampleList); do

        export ID
        CLEAN_UP_DEP=()

        module=0
        if [[ $count -ge 2 ]]; then continue; fi

        ##- 0.0 Copy raw Nanopore file
            module_setup 0.0_raw_reads.sh

            # Module inputs -------------
            header3="Copy raw Nanopore file"
            DEPENDENT_JOB=(COMPLETE)
            hpc_opts="--partition=Orion --nodes=1 --cpus-per-task=4 --mem=8GB --time=24:00:00 --mail-user=${email} --mail-type=FAIL"
            export raw_readDir
            # Skip copy if pre-QC is already done
            Complete_tag=(${pre_qcDir}/COMPLETE/${ID})
            # -----------------------------------

            run_module
            RAW_READS_JOB=$Current_Job


        ##- 0.1 Pre-QC (FastQC on raw reads)
            if $run_pre_qc; then
                module_setup 0.1_pre_qc.sh

                # Module inputs -------------
                header3="Pre-QC"
                DEPENDENT_JOB=(${RAW_READS_JOB##* })
                hpc_opts=$pre_qc_opts
                export raw_readDir
                export pre_qcDir
                # ---------------------------

                run_module
                PRE_QC_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
            else
                PRE_QC_JOB=COMPLETE
            fi


        # NOTE: 0.2 Deduplication is not run for ONT (ligation library prep has no PCR duplicates).


        ##- 0.3 Sequence trim (Porechop adapter trimming + NanoFilt quality/length filtering)
            if $run_trim; then
                module_setup 0.3_sequence_trim.sh

                # Module inputs -------------
                header3="Sequence trim (Porechop + NanoFilt)"
                DEPENDENT_JOB=(${PRE_QC_JOB##* })
                hpc_opts=$trim_opts
                export raw_readDir
                export trimmed_Dir
                # ---------------------------

                run_module
                TRIM_JOB=$Current_Job
            else
                TRIM_JOB=COMPLETE
            fi


        ##- 0.4 Host decontamination (minimap2 map-ont + samtools)
            if $run_decontam; then
                module_setup 0.4_host_decontamination.sh

                # Module inputs -------------
                header3="Host decontamination (minimap2 map-ont)"
                DEPENDENT_JOB=(${TRIM_JOB##* })
                hpc_opts=$decontam_opts
                export trimmed_Dir
                export clean_readDir
                export GRCh38_FASTA
                # ---------------------------

                run_module
                DECONTAM_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
            else
                DECONTAM_JOB=COMPLETE
            fi


        ##- 1.1 Assembly (metaFlye)
            if $run_asm; then
                module_setup 1.1_assembly.sh

                # Module inputs -------------
                header3="Assembly (metaFlye)"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$asm_opts
                export clean_readDir
                # ---------------------------

                run_module
                ASSEMBLY_JOB=$Current_Job
            else
                ASSEMBLY_JOB=COMPLETE
            fi


        ##- 1.2 Evaluation (filter contigs <${min_contig_len}bp)
            if $run_eval; then
                module_setup 1.2_evaluation.sh

                # Module inputs -------------
                header3="Evaluation (filter short contigs)"
                DEPENDENT_JOB=(${ASSEMBLY_JOB##* })
                hpc_opts=$eval_opts
                export assemblyDir
                export evaluationDir
                export min_contig_len
                # ---------------------------

                run_module
                EVALUATION_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
            else
                EVALUATION_JOB=COMPLETE
            fi


        ##- 2.1 Kraken2 (ONT reads + assembly)
            if $run_k2; then
                module_setup 2.1_kraken2.sh

                # Module inputs -------------
                header3="Kraken2 (ONT reads + assembly)"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$k2_opts
                export clean_readDir
                export evaluationDir
                if [[ ! -d 2.2_bracken ]]; then mkdir -p 2.2_bracken; fi
                # ---------------------------

                run_module
                KRAKEN_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
            else
                KRAKEN_JOB=COMPLETE
            fi


        ##- 3.1 Binning (MetaBat2 + MaxBin2 + CONCOCT, minimap2 -x map-ont alignment)
            if $run_bin; then
                module_setup 3.1_binning.sh

                # Module inputs -------------
                header3="Binning (MetaBat2 + MaxBin2 + CONCOCT)"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$bin_opts
                export clean_readDir
                export evaluationDir
                export binningDir
                # ---------------------------

                run_module
                BINNING_JOB=$Current_Job
            else
                BINNING_JOB=COMPLETE
            fi


        ##- 3.2 Refine bins (metaWRAP bin_refinement + CheckM)
            if $run_refine; then
                module_setup 3.2_refine_bins.sh

                # Module inputs -------------
                header3="Bin Refinement (metaWRAP + CheckM)"
                DEPENDENT_JOB=(${BINNING_JOB##* })
                hpc_opts=$refine_opts
                export binningDir
                export min_completion
                export max_contam
                # ---------------------------

                run_module
                REFINE_JOB=$Current_Job
            else
                REFINE_JOB=COMPLETE
            fi


        ##- 3.3 Reassemble bins (Flye)
            if $run_reassem; then
                module_setup 3.3_reassemble_bins.sh

                # Module inputs -------------
                header3="Bin Reassembly (Flye)"
                DEPENDENT_JOB=(${REFINE_JOB##* })
                hpc_opts=$reassem_opts
                export clean_readDir
                export refinedbinDir
                export min_completion
                export max_contam
                # ---------------------------

                run_module
                REASSEMBLY_JOB=$Current_Job
            else
                REASSEMBLY_JOB=COMPLETE
            fi


        ##- 4.1 Classify bins (GTDBtk)
            if $run_classify; then
                module_setup 4.1_classify_bins.sh

                # Module inputs -------------
                header3="Classify Bins (GTDBtk)"
                DEPENDENT_JOB=(${REASSEMBLY_JOB##* })
                hpc_opts=$classify_opts
                export reassemDir
                # ---------------------------

                run_module
                CLASSIFY_JOB=$Current_Job
            else
                CLASSIFY_JOB=COMPLETE
            fi


        ##- 4.2 Annotate bins (Prokka/Bakta)
            if $run_annotate; then
                module_setup 4.2_annotate_bins.sh

                # Module inputs -------------
                header3="Annotate Bins (Prokka/Bakta)"
                DEPENDENT_JOB=(${REASSEMBLY_JOB##* })
                hpc_opts=$annotate_opts
                export reassemDir
                # ---------------------------

                run_module
                ANNOTATE_JOB=$Current_Job
            else
                ANNOTATE_JOB=COMPLETE
            fi


        ##- 5.1 NT AMR from assembly (AMRFinder+ + RGI)
            if $run_amr_nt_asm; then
                module_setup 5.1_NT_amr_assembly.sh

                # Module inputs -------------
                header3="NT AMR from Assembly (AMRFinder+ + RGI)"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$amr_nt_opts
                export evaluationDir
                # ---------------------------

                run_module
                AMR_NT_ASM_JOB=$Current_Job
            else
                AMR_NT_ASM_JOB=COMPLETE
            fi


        ##- 5.4 RGI BWT (read-based AMR, single-end ONT)
            if $run_rgi_bwt; then
                module_setup 5.4_rgi_bwt.sh

                # Module inputs -------------
                header3="RGI BWT (read-based AMR, ONT single-end)"
                DEPENDENT_JOB=(${DECONTAM_JOB##* })
                hpc_opts=$rgi_bwt_opts
                export clean_readDir
                export rgi_bwt_dir
                # ---------------------------

                run_module
                RGI_BWT_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
            else
                RGI_BWT_JOB=COMPLETE
            fi


        ##- 5.5 AA AMR from assembly (Prodigal + Bakta + AMRFinder+ + RGI proteins)
            if $run_amr_aa_asm; then
                module_setup 5.5_AA_amr_assembly.sh

                # Module inputs -------------
                header3="AA AMR from Assembly (Prodigal + Bakta + AMRFinder+ + RGI)"
                DEPENDENT_JOB=(${EVALUATION_JOB##* })
                hpc_opts=$asm_profile_opts
                export evaluationDir
                export reassemDir
                # ---------------------------

                run_module
                AMR_AA_ASM_JOB=$Current_Job
                CLEAN_UP_DEP+=(${Current_Job##* })
            else
                AMR_AA_ASM_JOB=COMPLETE
            fi


        ##- COMPLETE (cleanup)
            if $run_clean_up; then
                module_setup COMPLETE.sh

                # Module inputs -------------
                header3="Pipeline cleanup"
                DEPENDENT_JOB=(${CLEAN_UP_DEP[@]})
                hpc_opts=$discard_opts
                # ---------------------------

                run_module
                COMPLETE_JOB=$Current_Job
            fi


    done
