#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Performs bin refinement using metaWRAP's bin_refinement module:
    #     1. Copies bins between 50kb and 20Mb.
    #     2. Performs bin refinement (binning_refiner.py).
    #     3. Runs CheckM on all sets of refined bins (checkm lineage_wf).
    #     4. Consolidates bin sets (consolidate_two_sets_of_bins.py).
    #     5. Dereplicates contigs between bins (dereplicate_contigs_in_bins.py).
    #     6. Runs CheckM on refined, consolidated, and dereplicated bins (checkm lineage_wf).
    #     7. Makes completion and contamination ranking plots for all refinement iterations (plot_binning_results.py).
    #     8. Makes completion and contamination ranking plots for the final output (plot_binning_results.py).

    # It is important to note that it has been designed for a specific working directory. 
    # Therefore, the reproduction of the results will require small modifications of the script 
    # or the adaptation of your working directory.

    # Created on Nov 7, 2024

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
    H1 "Description: 3.2_refine_bins.sh"
        echo -e "This script does the following:"
        echo -e "1. Performs bin refinement using metaWRAP's bin_refinement module:"
        echo -e "    1. Copies bins between 50kb and 20Mb."
        echo -e "    2. Performs bin refinement (binning_refiner.py)."
        echo -e "    3. Runs CheckM on all sets of refined bins (checkm lineage_wf)."
        echo -e "    4. Consolidates bin sets (consolidate_two_sets_of_bins.py)."
        echo -e "    5. Dereplicates contigs between bins (dereplicate_contigs_in_bins.py)."
        echo -e "    6. Runs CheckM on refined, consolidated, and dereplicated bins (checkm lineage_wf)."
        echo -e "    7. Makes completion and contamination ranking plots for all refinement iterations (plot_binning_results.py)."
        echo -e "    8. Makes completion and contamination ranking plots for the final output (plot_binning_results.py)."

    H1 "Job Context"
        OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
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

    H1 "Variables"
        comment "SampleID (ID): ${ID}"
        c=$min_completion; echo -e "Minimum completion: ${c}%"
        x=$max_contam; echo -e "Maximum contamination: ${x}%"

        H2 "Input"
            metabat2_OUT=${binningDir}/${ID}/metabat2_bins; echo -e "${metabat2_OUT}"
            maxbin2_OUT=${binningDir}/${ID}/maxbin2_bins; echo -e "${maxbin2_OUT}"
            concoct_OUT=${binningDir}/${ID}/concoct_bins; echo -e "${concoct_OUT}"

        H2 "Output"
            out=${moduleDir}/$ID
            if [[ ! -d ${out} ]]; then mkdir -p $out; fi

            echo -e "${out}/metawrap_${c}_${x}_bins"
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()


STEP="Bin Refinement"
    H1 "Bin refinement (>${c}% completeness, <${x}% contamination)"

    outputFile=${out}/metawrap_${c}_${x}_bins.stats
    if [[ -s "$outputFile" ]]; then
        comment "Output file already found. Skipping this command..."
    else
        start=$SECONDS
        #------------

        module load anaconda3/2023.09
        source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
        conda init
        conda activate metawrap-env

        metawrap bin_refinement \
            -o ${out} \
            -t $SLURM_CPUS_PER_TASK \
            -A ${binningDir}/${ID}/metabat2_bins \
            -B ${binningDir}/${ID}/maxbin2_bins \
            -C ${binningDir}/${ID}/concoct_bins \
            -c $c -x $x

        conda deactivate
        module unload anaconda3/2023.09

        #------------
        end=$SECONDS; duration=$(( end-start ))

        if [[ -s "$outputFile" ]]; then
            H2 "Yipee! $STEP Complete"
            comment "$STEP: $(elapsed_time "$duration")"
        fi
    fi
    step_completion "${outputFile}"


echo -e "\nmetaWRAP bins:"; awk -v com="$c" -v con="$x" '$2 > com && $3 < con' ${out}/metawrap_${c}_${x}_bins.stats | wc -l
echo -e "\nMetaBat2 bins:"; awk -v com="$c" -v con="$x" '$2 > com && $3 < con' ${out}/metabat2_bins.stats | wc -l
echo -e "\nMaxBin2 bins:"; awk -v com="$c" -v con="$x" '$2 > com && $3 < con' ${out}/maxbin2_bins.stats | wc -l
echo -e "\nCONCOCT bins:"; awk -v com="$c" -v con="$x" '$2 > com && $3 < con' ${out}/concoct_bins.stats | wc -l

# Completion
    module_completion

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
