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
    H1 "Description: 5.1_NT_amr_detection.sh"
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
        inputType=bin

        H2 "Input"
            if [[ $inputType == "asm" ]]; then
                inputFile=${evaluationDir}/${ID}_final_assembly.fasta
                echo -e $inputFile
            fi

            if [[ $inputType == "bin" ]]; then
                echo -e "${reassemDir}/${ID}/reassembled_bins"
            fi


        H2 "Output"
            if [[ ! -d ${moduleDir}/${ID} ]]; then mkdir -p ${moduleDir}/${ID} ; fi

            prodigal_output=${moduleDir}/${ID}/${ID}_genes.faa; echo -e "${prodigal_output}"
            bakta_output=${moduleDir}/${ID}/bakta/${ID}.tsv; echo -e "${bakta_output}"
            amr_output=${moduleDir}/${ID}/${ID}.amrfinder.tsv; echo -e "${amr_output}"
            rgi_output=${moduleDir}/$ID/${ID}.rgi.txt; echo -e "${rgi_output}"
            gene_list=${moduleDir}/$ID/${ID}_geneslist.tsv; echo -e "${gene_list}"
            outputFile=${moduleDir}/$ID/${ID}_gene_annotations.tsv; echo -e "${outputFile}"

            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()

# Load environments --------------------------------------------------------------------------------------------------------
    module load anaconda3/2023.09
    source $miniforge_init
    # source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    # conda activate metawrap-env



# Run functions ------------------------------------------------------------------------------------------------------------

func="Prodigal"
    H1 "$func"

    
    if [[ -s "$prodigal_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        if [[ $inputType == "asm" ]]; then
            prodigal -i ${evaluationDir}/${ID}_final_assembly.fasta -a ${moduleDir}/${ID}/${ID}_genes.faa -q
            if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi
        fi

        if [[ $inputType == "bin" ]]; then
            
            cat "${reassemDir}/${ID}/reassembled_bins"/*.fa > "${moduleDir}/${ID}/${ID}_all_reassembled_bins.fa"

            prodigal -i ${moduleDir}/${ID}/${ID}_all_reassembled_bins.fa -a ${moduleDir}/${ID}/${ID}_genes_original.faa -q

            in=${moduleDir}/${ID}/${ID}_genes_original.faa
            out=${moduleDir}/${ID}/${ID}_genes.faa
            map=${moduleDir}/${ID}/${ID}_genes.rename.map

            awk -v MAP="$map" '
              BEGIN{ OFS=""; }
              /^>/{
                full=$0
                name=substr($0,2); sub(/ .*/,"",name)          # header name up to first space
                if (++c[name] > 1){
                  new = name "_dup" c[name]
                  print name "\t" new >> MAP
                  print ">" new
                } else {
                  print name "\t" name >> MAP
                  print full
                }
                next
              }
              { print }
            ' "$in" > "$out"

            if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi
            substep_completion "$in" "$map"
        fi

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$prodigal_output" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi
    substep_completion "${prodigal_output}"


func="Bakta"
    H1 "$func"
    module load anaconda3/2023.09
    source $conda_init

    
    bakta_inf=${moduleDir}/${ID}/bakta/${ID}.inference.tsv
    if [[ -s "$bakta_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        if [[ -d ${moduleDir}/${ID}/bakta ]]; then rm -r ${moduleDir}/${ID}/bakta; fi

        conda activate $BAKTA_ENV
        export BAKTA_DB

        bakta_proteins --db ${BAKTA_DB} \
        --output ${moduleDir}/${ID}/bakta \
        --prefix ${ID} \
        --threads $OMP_NUM_THREADS \
        ${moduleDir}/${ID}/${ID}_genes.faa 

        if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi

        conda deactivate


        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$bakta_output" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${bakta_output}"


func="AMRFinderPlus"
    H1 "$func"
    source $miniforge_init

    # amr_output=${argDetect}/AMR/${ID}.amrfinder.txt
    
    if [[ -s "$amr_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load blast/2.11.0+
        module load hmmer/3.3.2
        # conda init
        conda activate $AMRFINDER_ENV

        amrfinder \
        --threads $SLURM_CPUS_PER_TASK \
        -p ${moduleDir}/${ID}/${ID}_genes.faa \
        -o ${moduleDir}/${ID}/${ID}.amrfinder.tsv

        if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi

        conda deactivate
        module unload blast/2.11.0+
        module unload hmmer/3.3.2

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$amr_output" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${amr_output}"


func="RGI"
    H1 "$func"
    # rgi_output=${argDetect}/RGI/${ID}.rgi.txt
    rgi_output=${moduleDir}/$ID/${ID}.rgi.txt
    if [[ -s "$rgi_output" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load diamond/2.0.9
        source $RGI_ENV

        tr -d "*" < ${moduleDir}/${ID}/${ID}_genes.faa > ${moduleDir}/${ID}/${ID}_genes_rgi.faa
        substep_completion "${moduleDir}/${ID}/${ID}_genes_rgi.faa"

        rgi main \
        -i ${moduleDir}/${ID}/${ID}_genes_rgi.faa \
        -t protein \
        -o ${moduleDir}/${ID}/${ID}.rgi \
        --clean \
        -n $SLURM_CPUS_PER_TASK

        if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi
        substep_completion "${moduleDir}/${ID}/${ID}.rgi.json"

        conda deactivate
        module unload diamond/2.0.9

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$rgi_output" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${rgi_output}"


func="Gene list"
    H1 "$func"
    gene_list=${moduleDir}/$ID/${ID}_geneslist.tsv
    if [[ -s "$gene_list" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        grep "^>" ${moduleDir}/${ID}/${ID}_genes.faa | awk '{print $1}' | sed 's/>//' > ${moduleDir}/$ID/${ID}_geneslist.tsv

        if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$gene_list" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${gene_list}"


func="Annotate genes using AMRFinder, RGI and Bakta results"
    H1 "$func"
    
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

        if [[ $? -ne 0 ]]; then error "Something went wrong with $func. Exiting"; fi

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$outputFile" ]]; then
            H2 "Yipee! $func Complete"
            comment "$func: $(elapsed_time "$duration")"
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
