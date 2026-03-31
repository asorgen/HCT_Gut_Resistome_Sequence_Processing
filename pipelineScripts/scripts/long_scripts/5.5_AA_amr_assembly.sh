#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. ORF prediction on the filtered assembly with Prodigal
    # 2. Functional annotation of predicted proteins with Bakta
    # 3. AMR annotation of predicted proteins with AMRFinderPlus
    # 4. AMR annotation of predicted proteins with RGI
    # 5. Gene list extraction
    # 6. Merges AMRFinder, RGI, and Bakta results into a unified gene annotations table

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

    # Input: ${evaluationDir}/${ID}_final_assembly.fasta (nucleotide FASTA from 1.2_evaluation.sh)
    # Output: ${moduleDir}/${ID}/${ID}_gene_annotations.tsv

# Slurm Resource Options ---------------------------------------------------------------------------------------------------

    # Job partition (--partition=<partition_names>; -p <partition_names>; SBATCH_PARTITION) | Options: Orion, Nebula, Pisces
    # Job name (--job-name=<name>; -J <name>; SBATCH_JOB_NAME)
    # Path to file storing text output. (--output=<filename_pattern>; -o <name>; SBATCH_OUTPUT)
    # Node count required for the job (--nodes=<count>; -N <count>)
    # Request that cpus per task (--cpus-per-task=<ncpus>)
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

# Print description script information to log ----------------------------------------------------------------------------------------
    /bin/date
    H1 "Description: 5.5_AA_amr_assembly.sh"
        echo -e "This script does the following:"
        echo -e "1. ORF prediction with Prodigal"
        echo -e "2. Functional annotation with Bakta"
        echo -e "3. AMR annotation with AMRFinderPlus"
        echo -e "4. AMR annotation with RGI"
        echo -e "5. Gene list extraction"
        echo -e "6. Merges results into unified gene annotations table"

# Print job context and resources to log ----------------------------------------------------------------------------------------
    H1 "Job Context"
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))
        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

        echo
        print "----- Resources Requested -----"
        print "Nodes:            $SLURM_NNODES"
        print "Cores / node:     $SLURM_CPUS_PER_TASK"
        print "Total memory:     $Total_Gb Gb"
        print "Wall-clock time:  $JobTime"
        print "-------------------------------"

# Print variables to log ----------------------------------------------------------------------------------------
    H1 "Variables"
        comment "SampleID (ID): ${ID}"
        inputType=asm

    H2 "Input"
        inputFile=${evaluationDir}/${ID}_final_assembly.fasta
        echo -e $inputFile

    H2 "Output"
        if [[ ! -d ${moduleDir}/${ID} ]]; then mkdir -p ${moduleDir}/${ID}; fi

        prodigal_output=${moduleDir}/${ID}/${ID}_genes.faa; echo -e "${prodigal_output}"
        bakta_output=${moduleDir}/${ID}/bakta/${ID}.tsv; echo -e "${bakta_output}"
        amr_output=${moduleDir}/${ID}/${ID}.amrfinder.tsv; echo -e "${amr_output}"
        rgi_output=${moduleDir}/$ID/${ID}.rgi.txt; echo -e "${rgi_output}"
        gene_list=${moduleDir}/$ID/${ID}_geneslist.tsv; echo -e "${gene_list}"
        outputFile=${moduleDir}/$ID/${ID}_gene_annotations.tsv; echo -e "${outputFile}"

        if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()


# Prodigal ORF prediction --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env
    
    STEP="Prodigal"
    H1 "$STEP"

    if [[ -s "$prodigal_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        prodigal -i ${inputFile} -a ${moduleDir}/${ID}/${ID}_genes.faa -q
        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))
    fi
    step_completion "${prodigal_output}"


# Bakta annotation --------------------------------------------------------------------------------------------------------
    STEP="Bakta"
    H1 "$STEP"
    module load anaconda3/2023.09
    source $conda_init

    bakta_inf=${moduleDir}/${ID}/bakta/${ID}.inference.tsv
    if [[ -s "$bakta_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        if [[ -d ${moduleDir}/${ID}/bakta ]]; then rm -r ${moduleDir}/${ID}/bakta; fi

        echo -e "Activate the Bakta conda environment"; conda activate $BAKTA_ENV
        echo -e "Export the Bakta db path: $BAKTA_DB"; export BAKTA_DB
        CMD="bakta_proteins --db ${BAKTA_DB} --output ${moduleDir}/${ID}/bakta --prefix ${ID} --threads $SLURM_CPUS_PER_TASK ${moduleDir}/${ID}/${ID}_genes.faa"
        echo $CMD
        $CMD

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        conda deactivate

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$bakta_output" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${bakta_output}"


# AMRFinderPlus annotation --------------------------------------------------------------------------------------------------------
    STEP="AMRFinderPlus"
    H1 "$STEP"
    source $miniforge_init

    if [[ -s "$amr_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load blast/2.11.0+
        module load hmmer/3.3.2
        conda activate $AMRFINDER_ENV

        amrfinder \
        --threads $SLURM_CPUS_PER_TASK \
        -p ${moduleDir}/${ID}/${ID}_genes.faa \
        -o ${moduleDir}/${ID}/${ID}.amrfinder.tsv

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        conda deactivate
        module unload blast/2.11.0+
        module unload hmmer/3.3.2

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$amr_output" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${amr_output}"


# RGI annotation --------------------------------------------------------------------------------------------------------
    STEP="RGI"
    H1 "$STEP"
    rgi_output=${moduleDir}/$ID/${ID}.rgi.txt
    if [[ -s "$rgi_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load blast/2.11.0+
        source $RGI_ENV

        rgi load --card_json $CARD_DB/card.json

        tr -d "*" < ${moduleDir}/${ID}/${ID}_genes.faa > ${moduleDir}/${ID}/${ID}_genes_rgi.faa

        rgi main \
        -i ${moduleDir}/${ID}/${ID}_genes_rgi.faa \
        -t protein \
        -o ${moduleDir}/${ID}/${ID}.rgi \
        --clean \
        -n $SLURM_CPUS_PER_TASK

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        conda deactivate
        module unload blast/2.11.0+

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$rgi_output" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${rgi_output}"


# Gene list extraction --------------------------------------------------------------------------------------------------------
    STEP="Gene list"
    H1 "$STEP"
    gene_list=${moduleDir}/$ID/${ID}_geneslist.tsv
    if [[ -s "$gene_list" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        grep "^>" ${moduleDir}/${ID}/${ID}_genes.faa | awk '{print $1}' | sed 's/>//' > ${moduleDir}/$ID/${ID}_geneslist.tsv

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$gene_list" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${gene_list}"


#  Merge annotations --------------------------------------------------------------------------------------------------------
    STEP="Annotate genes using AMRFinder, RGI and Bakta results"
    H1 "$STEP" 

    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        python ${ScriptPath}/helper_scripts/00_processgenes.py \
        --genes $gene_list \
        --amr $amr_output \
        --rgi $rgi_output \
        --bakta $bakta_output \
        --output $outputFile

        if [[ $? -ne 0 ]]; then error "Something went wrong with $STEP. Exiting"; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))
        if [[ -s "$outputFile" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${outputFile}"

# Completion
    output_exists=$(test_for_output "${Complete_tag[@]}")
    if $output_exists; then
        touch ${moduleDir}/COMPLETE/$ID

        # Remove intermediate files
        for int_file in "${Intermediate_files[@]}"; do
            echo -e "rm $int_file"; rm $int_file
        done
    else
        error "Not all outputs were created."
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
