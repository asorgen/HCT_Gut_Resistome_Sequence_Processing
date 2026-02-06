#!/bin/bash
# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Performs functional annotation of bin sets with Prokka using metaWRAP's annotate_bins module:
    #     1. Shortens contig names to run Prokka (shorten_contig_names.py).
    #     2. Runs Prokka. 
    #     3. Pulls out functional annotations for each bin.
    #     4. Filters the translated genes.
    #     5. Filters the untranslated genes.

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 10, 2025

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
    H1 "Description: 5.4_rgi_bwt.sh"
        echo -e "This script does the following:"
        echo -e "1. Performs AMR gene annotation of .fna files using AMRFinderPlus"
        echo -e "2. Performs AMR gene annotation of .fna files using RGI"

    H1 "Job Context"
        OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
        print "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        print "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))

        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

        echo 
        print "----- Resources Requested -----"
        print "Nodes:            $SLURM_NNODES"
        print "Cores / node:     $SLURM_CPUS_PER_TASK"
        print "Total memory:     $Total_Gb Gb"
        print "Wall-clock time:  $JobTime"
        print "-------------------------------"

    H1 "Variables"
        comment "SampleID: ${ID}"

        H2 "Input"
            datasetROOT=$(pwd); #echo "$datasetROOT"
            READ=${datasetROOT}/${clean_readDir}/${ID}
            R1=${READ}_1.fastq.gz; echo $R1
            R2=${READ}_2.fastq.gz; echo $R2
            aligner=kma

        H2 "Output"
            if [[ ! -d ${moduleDir}/${aligner}_output ]]; then mkdir ${moduleDir}/${aligner}_output; fi
            outputFile=${moduleDir}/${aligner}_output/${ID}.rgi_${aligner}.txt
            echo -e "${outputFile}"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi

    # H2 "[ Start ]"
    # /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Run RGI-BWT --------------------------------------------------------------------------------------------------------------
    H1 "Running RGI-BWT"
    if [[ ! -f "${outputFile}" ]]; then # if the final output file doesn't exist 
        # Load environments
            comment "Load environments"
            module load anaconda3/2023.09
            module load diamond/2.0.9
            source $miniforge_init
            conda activate $RGI_ENV

        # Set up temp directory in scratch
            comment "$ID processing..."
            cd $CARD_DB
            mkdir -p ${ID}_${aligner}
            cd ${ID}_${aligner}

        # Load local CARD database
            if [[ ! -d "localDB" ]]; then # if the local database does not exist
                comment "Load CARD database"
                rgi load --card_json $CARD_DB/card.json --local
                rgi card_annotation -i $CARD_DB/card.json > card_annotation.log 2>&1
                rgi load -i $CARD_DB/card.json --card_annotation card_database_v4.0.1.fasta --local
            fi

        # Run rgi bwt
            if [[ ! -f "${ID}.gene_mapping_data.txt" ]]; then # if the gene mapping output does not exist
                comment "Run RGI-BWT for the first time."
                rgi bwt \
                --read_one ${R1} --read_two ${R2} \
                -a $aligner \
                --output_file ${ID} \
                --threads $SLURM_CPUS_PER_TASK \
                --debug \
                --clean --local &> ${ID}_log.out
                rgi_exit=$?
                if [[ $rgi_exit -ne 0 ]]; then
                    error "RGI-BWT failed with exit code $rgi_exit. Check ${ID}_log.out for details."
                fi

                if [[ $(cat ${ID}.gene_mapping_data.txt | wc -l) -gt 1 ]]; then
                    hits=$(cat ${ID}.gene_mapping_data.txt | wc -l); hits=$(( hits - 1 ))
                    comment "RGI-BWT has successfully completed."
                    comment "$hits hits"
                else
                    error "No hits found"
                fi
            elif [[ $(cat ${ID}.gene_mapping_data.txt | wc -l) -eq 1 ]]; then # else, if the gene mapping output does exist, but has no hits
                comment "Rerun RGI-BWT since no hits were found."
                rgi bwt \
                --read_one ${R1} --read_two ${R2} \
                -a $aligner \
                --output_file ${ID} \
                --threads $SLURM_CPUS_PER_TASK \
                --debug \
                --clean --local &> ${ID}_log.out
                rgi_exit=$?
                if [[ $rgi_exit -ne 0 ]]; then
                    error "RGI-BWT failed with exit code $rgi_exit. Check ${ID}_log.out for details."
                fi

                if [[ $(cat ${ID}.gene_mapping_data.txt | wc -l) -gt 1 ]]; then
                    hits=$(cat ${ID}.gene_mapping_data.txt | wc -l); hits=$(( hits - 1 ))
                    comment "RGI-BWT has successfully completed."
                    comment "$hits hits"
                else
                    error "No hits found"
                fi
            else
                comment "RGI-BWT has already successfully completed."
                hits=$(cat ${ID}.gene_mapping_data.txt | wc -l); hits=$(( hits - 1 ))
                comment "$hits hits"
            fi

        # Copy the gene mapping output to the project directory
            mkdir -p $datasetROOT/5.4_rgi_bwt/kma_output
            cp ${ID}.gene_mapping_data.txt $datasetROOT/5.4_rgi_bwt/kma_output/${ID}.rgi_kma.txt
            cd $datasetROOT

    else
        comment "RGI-BWT output has already been generated."
    fi

# Clean up temp output from scratch directory
    if [[ ! -f "$outputFile" ]] || [[ $(cat $outputFile | wc -l) -le 1 ]]; then
        if [[ -f "$outputFile" ]]; then rm $outputFile; fi
        error "RGI-BWT found no hits for this sample."

    else
        hits=$(cat $datasetROOT/5.4_rgi_bwt/kma_output/${ID}.rgi_kma.txt | wc -l); hits=$(( hits - 1 ))
        comment "$hits hits"
        touch ${moduleDir}/COMPLETE/$ID
        if [[ -d ${CARD_DB}/${ID}_${aligner} ]]; then 
            comment "Removing tmp scratch directory"
            rm -r ${CARD_DB}/${ID}_${aligner}
        fi
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
