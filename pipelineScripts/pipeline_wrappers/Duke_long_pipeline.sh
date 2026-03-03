#!/bin/bash

# Description
    #-------------------------------------------------------------------------------------------------------------------
    # This script is used to begin the HCT ARG Assembly pipeline for the Oxford Nanopore (long) read samples.

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
read=long
dataset=${cohort}_${read}

source ${HOME}/.bashrc
export pipelineConfig=${dataset}-read.config
source pipelineScripts/configs/${pipelineConfig}

if [[ ! -d $datasetDir ]]; then mkdir -p $datasetDir; fi
cd $datasetDir
mkdir -p LOGs

export seqPath
export readType

# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; }
    error () { echo $1; exit 1; }
    job_lookup() { squeue -u ${USER} --format='%.18i   %.9P   %.40j   %.1T   %.12M   %.10l   %.6D   %R' | awk -v id="$jobID" 'match($3,id) {print $1}'; }
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
        comment "[ Raw sequence directory ]: ${seqPath}"
        
fi


# Initialize count variable
count=0

export totalSamples=`tail -n +2 $sampleList | wc -l`

for ID in $(tail -n +2 $sampleList); do

    module=0
    export ID

    ##- Copy raw nanopore file
        
        module_setup 0.0_raw_reads.sh
        
        # Module specific -------------------
        header3="Copy raw nanopore file"
        DEPENDENT_JOB=COMPLETE
        hpc_opts=$copy_opts
        pipeline_tag=copy_ont
        # -----------------------------------
        
        # Module specific actions -----------
        Complete_tag=${moduleDir}/${ID}.fastq.gz
        Complete_tag=0.1_read_qc/${ID}/COMPLETE
        # -----------------------------------
        
        run_module 
        RAW_READS_JOB=$Current_Job       
        if [[ "$1" = "$pipeline_tag" ]]; then continue; fi
        


    ##- Read QC

        module_setup 0.1_read_qc.sh

        # Module inputs -------------
        header3="Read QC"
        DEPENDENT_JOB=COMPLETE
        hpc_opts=$read_qc_opts
        pipeline_tag=read_qc
        # ---------------------------


        # Module specific -----------
        export raw_readDir=0.0_raw_reads
        if [[ ! -d ${raw_readDir} ]]; then mkdir -p $raw_readDir; fi
        export clean_readDir=0.2_clean_reads
        if [[ ! -d ${clean_readDir} ]]; then mkdir -p $clean_readDir; fi
        mkdir -p ${moduleDir}/pre-QC_reports
        mkdir -p ${moduleDir}/post-QC_reports
        # ---------------------------
        
        run_module
        READ_QC_JOB=$Current_Job        
        if [[ "$1" = "$pipeline_tag" ]]; then continue; fi


    ##- Assembly
        
        module_setup 1.1_assembly.sh
        
        # Module specific -------------------
        header3="Assembly"
        DEPENDENT_JOB=$READ_QC_JOB
        hpc_opts=$asm_opts
        pipeline_tag=assembly
        # -----------------------------------
        
        # Module specific actions -----------
        # -----------------------------------
        
        
        run_module 
        ASSEMBLY_JOB=$Current_Job       
        if [[ "$1" = "$pipeline_tag" ]]; then continue; fi



    ##- Evaluation
        
        module_setup 1.2_evaluation.sh
        
        # Module specific -------------------
        header3="Evaluation"
        DEPENDENT_JOB=$ASSEMBLY_JOB
        hpc_opts=$eval_opts
        pipeline_tag=evaluation
        # -----------------------------------
        
        # Module specific actions -----------
        Complete_tag=${moduleDir}/${ID}_final_assembly.fasta
        # -----------------------------------
        
        run_module 
        EVALUATION_JOB=$Current_Job       
        if [[ "$1" = "$pipeline_tag" ]]; then continue; fi



    ##- Kraken2
        
        module_setup 2.1_kraken2.sh
        
        # Module specific -------------------
        header3="Kraken2"
        DEPENDENT_JOB=$EVALUATION_JOB
        hpc_opts=$k2_opts
        pipeline_tag=kraken2
        # -----------------------------------
        
        # Module specific actions -----------
        if [[ ! -d 2.2_bracken ]]; then mkdir -p 2.2_bracken; fi
        # -----------------------------------
        
        run_module 
        KRAKEN_JOB=$Current_Job       
        if [[ "$1" = "$pipeline_tag" ]]; then continue; fi

    # ##- Binning
        
    #     script=binning.sh
    #     moduleDir="${script%.sh}"
    #     export binningDir=$moduleDir

    #     DEPENDENT_JOB=$EVALUATION_JOB

    #     binning_OUT=${moduleDir}/${ID}/COMPLETE

    #     if [[ -a "$binning_OUT" ]]; then
    #         BINNING_JOB=COMPLETE
    #     else
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then

    #             if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
    #                 count=$((count + 1))
    #                 H2 "${ID}"
    #                 H3 "Binning"
    #                 BINNING_JOB=$(sbatch \
    #                     $bin_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                 echo $BINNING_JOB
    #                 echo "[ Log file ] -> LONG/${logDir}/${ID}.${BINNING_JOB##* }.log"
    #             else
    #                 H3 "Binning"
    #                 BINNING_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
    #                     $bin_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                 echo $BINNING_JOB
    #                 echo "[ Log file ] -> LONG/${logDir}/${ID}.${BINNING_JOB##* }.log"
    #             fi
    #         else
    #             BINNING_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "binning" ]]; then continue; fi

    # ##- Bin Refinement

    #     script=refine_bins.sh
    #     moduleDir="${script%.sh}"
    #     export refinedDir=$moduleDir

    #     export max_completion
    #     export min_contam

        
    #     refine_OUT=${moduleDir}/${ID}/COMPLETE

    #     if [[ -a "$refine_OUT" ]]; then
    #         REFINE_JOB=COMPLETE
    #     else
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$BINNING_JOB" = "COMPLETE" ]]; then
    #                 count=$((count + 1))
    #                 H2 "${ID}"
    #                 H3 "Bin Refinement"
    #                 REFINE_JOB=$(sbatch \
    #                     $refine_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $REFINE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${REFINE_JOB##* }.log"
    #             else
    #                 H3 "Bin Refinement"
    #                 REFINE_JOB=$(sbatch --dependency=afterok:${BINNING_JOB##* } \
    #                     $refine_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $REFINE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${REFINE_JOB##* }.log"
    #            fi
    #         else
    #             REFINE_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "refine" ]]; then continue; fi

    # ##- Bin Reassembly

    #     script=reassemble_bins.sh
    #     moduleDir="${script%.sh}"
    #     export reassemDir=$moduleDir

    #     DEPENDENT_JOB=$REFINE_JOB
    #     reassembly_OUT=${moduleDir}/${ID}/COMPLETE

    #     if [[ -a "$reassembly_OUT" ]]; then
    #         REASSEMBLY_JOB=COMPLETE
    #     else
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
    #                 count=$((count + 1))
    #                 H2 "${ID}"
    #                 H3 "Bin Reassembly"
    #                 REASSEMBLY_JOB=$(sbatch \
    #                     $reassem_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $REASSEMBLY_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${REASSEMBLY_JOB##* }.log"
    #             else
    #                 H3 "Bin Reassembly"
    #                 REASSEMBLY_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
    #                     $reassem_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $REASSEMBLY_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${REASSEMBLY_JOB##* }.log"
    #            fi
    #         else
    #             REASSEMBLY_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "reassembly" ]]; then continue; fi

    # ##- Bin Classification

    #     script=classify_bins.sh
    #     moduleDir="${script%.sh}"

    #     DEPENDENT_JOB=$REASSEMBLY_JOB

    #     classify_OUT=${moduleDir}/${ID}/COMPLETE

    #     if [[ -a "$classify_OUT" ]]; then
    #         CLASSIFY_JOB=COMPLETE
    #     else
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
    #                 count=$((count + 1))
    #                 H2 "${ID}"
    #                 H3 "Bin Classification"
    #                 CLASSIFY_JOB=$(sbatch \
    #                     $classify_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $CLASSIFY_JOB
    #                echo "[ Log file ] -> ${logDir}/${ID}.${CLASSIFY_JOB##* }.log"
    #             else
    #                 H3 "Bin Classification"
    #                 CLASSIFY_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
    #                     $classify_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $CLASSIFY_JOB
    #                echo "[ Log file ] -> ${logDir}/${ID}.${CLASSIFY_JOB##* }.log"
    #            fi
    #         else
    #             CLASSIFY_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "classify" ]]; then continue; fi

    # ##- Functional Annotation

    #     script=annotate_bins.sh
    #     moduleDir="${script%.sh}"

    #     DEPENDENT_JOB=$REASSEMBLY_JOB
    #     annotate_OUT=${moduleDir}/${ID}/COMPLETE

    #     if [[ -f "$annotate_OUT" ]]; then
    #         ANNOTATE_JOB=COMPLETE
    #     else
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
    #                 count=$((count + 1))
    #                 H2 "${ID}"
    #                 H3 "Functional Annotation"
    #                 ANNOTATE_JOB=$(sbatch \
    #                     $annotate_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $ANNOTATE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${ANNOTATE_JOB##* }.log"
    #             else
    #                 H3 "Functional Annotation"
    #                 ANNOTATE_JOB=$(sbatch --dependency=afterok:${DEPENDENT_JOB##* } \
    #                     $annotate_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $ANNOTATE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${ANNOTATE_JOB##* }.log"
    #            fi
    #         else
    #             ANNOTATE_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "annotation" ]]; then continue; fi

    # ##- AMR Identification

    #     script=amr_detection.sh
    #     moduleDir="${script%.sh}"

    #     mkdir -p ${moduleDir}/AMR
    #     mkdir -p ${moduleDir}/RGI
    #     amr_OUT=${moduleDir}/RGI/${ID}.rgi.txt

    #     if [[ -s "$amr_OUT" ]]; then
    #         AMR_JOB=COMPLETE
    #     else
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$EVALUATION_JOB" = "COMPLETE" ]]; then
    #                 count=$((count + 1))
    #                 H2 "${ID}"
    #                 H3 "AMR Detection"
    #                 AMR_JOB=$(sbatch \
    #                     $amr_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $AMR_JOB
    #                echo "[ Log file ] -> ${logDir}/${ID}.${AMR_JOB##* }.log"
    #             else
    #                 H3 "AMR Detection"
    #                 AMR_JOB=$(sbatch --dependency=afterok:${EVALUATION_JOB##* } \
    #                     $amr_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                 echo $AMR_JOB
    #                 echo "[ Log file ] -> ${logDir}/${ID}.${AMR_JOB##* }.log"
    #            fi
    #         else
    #             AMR_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "amr_search" ]]; then continue; fi
    
    # ##- Assembly Gene Profiling
    
    #     inputType=asm
    #     script=gene_profiling.sh
    #     moduleDir="${inputType}_${script%.sh}"
    #     DEPENDENT_JOB=$EVALUATION_JOB

    #     mkdir -p ${moduleDir}/
    #     asm_OUT=${moduleDir}/$ID/${ID}_gene_annotations.tsv

    #     if [[ -s "$asm_OUT" ]]; then
    #         ASM_PROFILE_JOB=COMPLETE
    #     else
            
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export inputType
    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
    #                 if [[ "$AMR_JOB" = "COMPLETE" ]]; then count=$((count + 1)); H2 "${ID}"; fi
                    
    #                 H3 "Assembly Gene Profiling"
    #                 ASM_PROFILE_JOB=$(sbatch \
    #                     $asm_profile_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $ASM_PROFILE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${ASM_PROFILE_JOB##* }.log"
    #             else
    #                 H3 "Assembly Gene Profiling"
    #                 ASM_PROFILE_JOB=$(sbatch --dependency=afterok:${EVALUATION_JOB##* } \
    #                     $asm_profile_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $ASM_PROFILE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${ASM_PROFILE_JOB##* }.log"
    #            fi
    #         else
    #             ASM_PROFILE_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "asm_profile" ]]; then continue; fi
    
    # ##- Bin Gene Profiling
    
    #     inputType=bin
    #     script=gene_profiling.sh
    #     moduleDir="${inputType}_${script%.sh}"

    #     DEPENDENT_JOB=$EVALUATION_JOB

    #     mkdir -p ${moduleDir}/
    #     bin_OUT=${moduleDir}/$ID/${ID}_gene_annotations.tsv

    #     if [[ -s "$bin_OUT" ]]; then
    #         BIN_PROFILE_JOB=COMPLETE
    #     else
            
    #         jobID=$(printf "%02d" $module)
    #         jobID=${jobID}_${ID}_${readType}

    #         export inputType
    #         export moduleDir
    #         logDir=${moduleDir}/logs
    #         mkdir -p $logDir

    #         queued=$(job_lookup $jobID)
    #         if [ -z "$queued" ]; then
    #             if [[ "$DEPENDENT_JOB" = "COMPLETE" ]]; then
    #                 if [[ "$ASM_PROFILE_JOB" = "COMPLETE" ]]; then count=$((count + 1)); H2 "${ID}"; fi
    #                 H3 "Bin Gene Profiling"
    #                 BIN_PROFILE_JOB=$(sbatch \
    #                     $asm_profile_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $BIN_PROFILE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${BIN_PROFILE_JOB##* }.log"
    #             else
    #                 H3 "Bin Gene Profiling"
    #                 BIN_PROFILE_JOB=$(sbatch --dependency=afterok:${EVALUATION_JOB##* } \
    #                     $asm_profile_opts \
    #                     --job-name=${jobID} -o ${logDir}/${ID}.%A.log ./scripts/modules/${script})
    #                echo $BIN_PROFILE_JOB
    #                echo "[ Log file ] -> LONG/${logDir}/${ID}.${BIN_PROFILE_JOB##* }.log"
    #            fi
    #         else
    #             BIN_PROFILE_JOB=$queued
    #         fi
    #     fi
    #     ((module++))

    #     if [[ "$1" = "bin_profile" ]]; then continue; fi

    # # if [[ $count -ge 1 ]]; then exit 0; fi

done


