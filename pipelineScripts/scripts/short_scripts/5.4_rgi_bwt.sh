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

# Set function for output comments -----------------------------------------------------------------------------------------
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }

# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 5.4_rgi_bwt.sh"
        echo -e "This script does the following:"
        echo -e "1. Performs AMR gene annotation of .fna files using AMRFinderPlus"
        echo -e "2. Performs AMR gene annotation of .fna files using RGI"

    H1 "Job Context"
        OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
        comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
        comment "Running on host: `hostname`"

        Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))

        JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

        echo 
        comment "----- Resources Requested -----"
        comment "Nodes:            $SLURM_NNODES"
        comment "Cores / node:     $SLURM_CPUS_PER_TASK"
        comment "Total memory:     $Total_Gb Gb"
        comment "Wall-clock time:  $JobTime"
        comment "-------------------------------"

    H1 "Variables"
        comment "SampleID (ID): ${ID}"

        H2 "Input"
            datasetROOT=$(pwd); #echo "$datasetROOT"
            READ=${datasetROOT}/${clean_readDir}/${ID}
            R1=${READ}_1.fastq.gz; echo $R1
            R2=${READ}_2.fastq.gz; echo $R2
            aligner=kma


        H2 "Output"
            outputFile=${moduleDir}/${aligner}_output/${ID}.rgi_${aligner}.txt
            echo -e "${outputFile}"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()




if [[ ! -f "${outputFile}" ]]; then 

    # miniforge_init=/users/asorgen/miniforge3/etc/profile.d/conda.sh
    # RGI_ENV=/users/asorgen/miniforge3/envs/rgi
    # CARD_DB=/scratch/asorgen/RGI_databases/rgi_4.0.1
    datasetROOT=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Pipeline/sequence_processing/UNC_short
    datasetROOT=${ROOT}/${dataset}
    # clean_readDir=0.4_host_decontamination

    cd $CARD_DB

    module load anaconda3/2023.09
    module load diamond/2.0.9
    source $miniforge_init
    conda activate $RGI_ENV

    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }

    # ID=BMT101D-7


    READ=${datasetROOT}/${clean_readDir}/${ID}
    reads_1=${READ}_1.fastq.gz
    reads_2=${READ}_2.fastq.gz
    aligner=kma

    comment "$ID processing..."
    mkdir -p ${ID}_${aligner}
    cd ${ID}_${aligner}

    # trimmed_reads_1=/scratch/asorgen/RGI_databases/rgi_4.0.1/${ID}_kma/trimgalore_output/${ID}_1_val_1.fq.gz
    # trimmed_reads_2=/scratch/asorgen/RGI_databases/rgi_4.0.1/${ID}_kma/trimgalore_output/${ID}_2_val_2.fq.gz

    # if [[ ! -e "${trimmed_reads_1}" || ! -e "${trimmed_reads_2}" ]]; then
    #     comment "$STEP"
    #     conda activate metawrap-env
    #     mkdir trimgalore_output
    #     trim_galore --paired $reads_1 $reads_2 --gzip --quality 30 --length 50 --stringency 3 --max_n 0 -o trimgalore_output/
    #     conda deactivate
    # fi


    conda activate rgi

    rgi load --card_json $CARD_DB/card.json --local
    rgi card_annotation -i $CARD_DB/card.json > card_annotation.log 2>&1
    rgi load -i $CARD_DB/card.json --card_annotation card_database_v4.0.1.fasta --local

    rgi bwt \
    --read_one ${reads_1} --read_two ${reads_2} \
    -a $aligner \
    --output_file ${ID} \
    --threads $SLURM_CPUS_PER_TASK \
    --debug \
    --clean --local

    cat ${ID}.gene_mapping_data.txt | wc -l

    cp ${ID}.gene_mapping_data.txt $datasetROOT/5.4_rgi_bwt/kma_output/${ID}.rgi_kma.txt
    cat $datasetROOT/5.4_rgi_bwt/kma_output/${ID}.rgi_kma.txt | wc -l

    cd ..
    rm -r ${ID}_kma
fi

if [[ $(cat $outputFile | wc -l) -eq 1 ]]; then 
    echo -e "rm $outputFile"
    rm $outputFile
    error "RGI BWT found no hits for this sample."

else
    touch ${moduleDir}/COMPLETE/$ID
fi


exit 0
# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    # source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    # conda activate metawrap-env


# Run functions ------------------------------------------------------------------------------------------------------------
if [[ $(cat $outputFile | wc -l) -eq 1 ]]; then 
    rm $outputFile
fi


func="RGI"
    H1 "$func"
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        STEP="Unzip the sequence files"
            run_step=false
            if $run_step; then
                reads_1=${moduleDir}/${ID}_1.fastq
                reads_2=${moduleDir}/${ID}_2.fastq
                Unzip_files=(${reads_1} ${reads_2})
                output_exists=$(test_for_output "${Unzip_files[@]}")
                if ! "$output_exists"; then
                    comment "$STEP"
                    gunzip -c $R1 > ${reads_1}
                    gunzip -c $R2 > ${reads_2}
                fi
                substep_completion "${Unzip_files[@]}"
            else
                reads_1=$R1
                reads_2=$R2
            fi
            run_step=true
     

        comment "Make a tmp directory for processing."
            comment "$ID processing..."
            mkdir -p $CARD_DB/${ID}_${aligner}
            cd $CARD_DB/${ID}_${aligner}

        STEP="Additional quality filtering with Trim Galore"

            trimmed_reads_1=$RGI_local/${ID}_${aligner}/trimgalore_output/${ID}_1_val_1.fq.gz
            trimmed_reads_2=$RGI_local/${ID}_${aligner}/trimgalore_output/${ID}_2_val_2.fq.gz

            if [[ ! -e "${trimmed_reads_1}" || ! -e "${trimmed_reads_2}" ]]; then
                comment "$STEP"
                conda activate metawrap-env
                mkdir trimgalore_output
                trim_galore --paired $reads_1 $reads_2 --gzip --quality 30 --length 50 --stringency 3 --max_n 0 -o trimgalore_output/
                conda deactivate
            fi


        comment "Load tools"
            module load diamond/2.0.9
            source $miniforge_init
            conda activate $RGI_ENV


        STEP="Load CARD reference data into local or working directory"
            CARD_version=$(jq '._version' $CARD_DB/card.json  | tr -d '"')
            H2 "CARD Database version: $CARD_version"
            # run_step=false
            if $run_step; then
            # if [[ ! -e "${RGI_local}/localDB/card_reference.fasta" ]]; then
                comment "$STEP"
                rgi load --card_json $CARD_DB/card.json --local
                       
                rgi card_annotation -i $CARD_DB/card.json > $CARD_DB/card_annotation.log 2>&1

                rgi load -i $CARD_DB/card.json --card_annotation $CARD_DB/card_database_v${CARD_version}.fasta
            fi

        STEP="Load CARD wildcard reference data into local or working directory"
            run_step=false
            if $run_step; then
            # if [[ ! -e "${RGI_local}/localDB/card_wildcard_reference.fasta" ]]; then
                comment "$STEP"
                rgi wildcard_annotation -i wildcard_v${CARD_version}--card_json card.json \
                    -v ${CARD_version} > wildcard_annotation.log 2>&1
                rgi load --wildcard_annotation wildcard_database_v${CARD_version}.fasta \
                    --card_json card.json \
                    --wildcard_index wildcard_v${CARD_version}index-for-model-sequences.txt \
                    --card_annotation card_database_v${CARD_version}.fasta --local
            fi

        STEP="Test kma manually"
            run_step=false
            if $run_step; then
                H2 "$STEP"
                kma \
                    -ipe \
                    ${reads_1} \
                    ${reads_2} \
                    -t_db ${RGI_local}/localDB/bwt/card_reference/kma/card_reference\
                    -o ${ID}_kma_test \
                    -t $SLURM_CPUS_PER_TASK
            fi
            run_step=true

        STEP="Downsample reads for time reference."
            run_step=false
            if $run_step; then
                comment "$STEP"
                echo $(( $(zcat "$reads_1" | wc -l) / 4 ))
                echo $(( $(zcat "$reads_2" | wc -l) / 4 ))

                zcat ${reads_1} | head -n 1000000 > ${ID}_R1_1M_subsample.fastq
                zcat ${reads_2} | head -n 1000000 > ${ID}_R2_1M_subsample.fastq
                reads_1=${ID}_R1_1M_subsample.fastq
                reads_2=${ID}_R2_1M_subsample.fastq
            fi
            run_step=true
        
        STEP="Run rgi bwt on ${READ}_1.fastq.gz and ${READ}_2.fastq.gz"
            # run_step=false
            if $run_step; then
                H2 "$STEP"
                rgi bwt \
                    --read_one ${trimmed_reads_1} --read_two ${trimmed_reads_2} \
                    -a $aligner \
                    --output_file ${ID} \
                    --threads $SLURM_CPUS_PER_TASK \
                    --debug \
                    --clean --local
                if [[ ! -s ${RGI_local}/${ID}_${aligner}/$ID.gene_mapping_data.txt ]]; then error "${RGI_local}/${ID}_${aligner}/$ID.gene_mapping_data.txt not found. Exiting"; fi
            fi
            run_step=true

            cd $datasetROOT

        STEP="Move output file to $moduleDir"
            # run_step=false
            if $run_step; then
                comment "$STEP"
                cp ${RGI_local}/${ID}_${aligner}/$ID.gene_mapping_data.txt $moduleDir/${ID}.rgi_${aligner}.txt
            fi
            run_step=true


        STEP="Remove tmp files"
            run_step=false
            if $run_step; then
                comment "$STEP"
                rm -r ${RGI_local}/${ID}_${aligner}
                # rm ${RGI_local}/$ID*
            fi
            run_step=true
        
        
        conda deactivate
        module unload diamond/2.0.9

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$outputFile" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        else
            error "Something went wrong with $func. Exiting"
        fi
    fi

if [[ $(cat $outputFile | wc -l) -eq 1 ]]; then 
    echo -e "rm $outputFile"
    rm $outputFile
    error "RGI BWT found no hits for this sample."

else
    touch ${moduleDir}/COMPLETE/$ID
fi

H1 "MODULE COMPLETE"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
