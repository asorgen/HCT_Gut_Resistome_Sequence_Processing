#!/bin/bash

# Description --------------------------------------------------------------------------------------------------------------
    # This script does the following:
    # 1. Recruits reads to bins for reassembly (minimap2).
    # 2. Filters recruited reads.
    # 3. Reassembles bins using Flye.
    # 4. Removes short contigs from reassemblies.
    # 5. Runs CheckM on all reassembled and original bins.
    # 6. Selects the best version of each bin (reassembled or original).
    # 7. Re-runs CheckM on the final best bins.
    # 8. Generates QA plots and summary stats.

    # It is important to note that it has been designed for a specific working directory.
    # Therefore, the reproduction of the results will require small modifications of the script
    # or the adaptation of your working directory.

    # Created on March 2026

    # @author: Alicia Sorgen - UNC Charlotte Dept of Bioinformatics and Genomics

    # Version: 1

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


# Print script information to log ------------------------------------------------------------------------------------------
    H1 "Description: 3.3_reassemble_bins.sh"
        echo -e "This script does the following:"
        echo -e "1. Recruits reads to bins for reassembly (minimap2)."
        echo -e "2. Filters recruited reads."
        echo -e "3. Reassembles bins using Flye."
        echo -e "4. Removes short contigs from reassemblies."
        echo -e "5. Runs CheckM on all reassembled and original bins."
        echo -e "6. Selects the best version of each bin (reassembled or original)."
        echo -e "7. Re-runs CheckM on the final best bins."
        echo -e "8. Generates QA plots and summary stats."

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

        H2 "Input"
            nanopore_reads=${clean_readDir}/${ID}_ont.fastq; echo -e "${nanopore_reads}"
            bins=${refinedbinDir}/${ID}/metawrap_${min_completion}_${max_contam}_bins; echo -e "${bins}"

        H2 "Output"
            out=${moduleDir}/${ID}
            if [[ ! -d ${out} ]]; then mkdir -p $out; fi
            echo -e "${out}/reassembled_bins.stats"
            if [[ ! -d ${moduleDir}/COMPLETE ]]; then mkdir -p ${moduleDir}/COMPLETE; fi

    H2 "[ Start ]"
    /bin/date
    SECONDS=0
    Complete_tag=()
    Intermediate_files=()


# metaWRAP helper functions ------------------------------------------------------------------------------------------------
    comm () { ${SOFT}/print_comment.py "$1" "-"; }
    warning () { ${SOFT}/print_comment.py "$1" "*"; }
    announcement () { ${SOFT}/print_comment.py "$1" "#"; }

# Parallelization functions ------------------------------------------------------------------------------------------------
    open_sem(){
        mkfifo pipe-$$
        exec 3<>pipe-$$
        rm pipe-$$
        local i=$1
        for((;i>0;i--)); do
            printf %s 000 >&3
        done
    }
    run_with_lock(){
        local x
        read -u 3 -n 3 x && ((0==x)) || exit $x
        (
        "$@"
        printf '%.3d' $? >&3
        )&
    }

# Default params -----------------------------------------------------------------------------------------------------------
    threads=1; mem=40; comp=70; cont=10; len=500
    bins=None; f_reads=None; r_reads=None; out=None

    # long options defaults
    strict_max=2; permissive_max=5
    run_checkm=true
    run_parallel=false
    nanopore=false
    mdmcleaner=false

# Load in params -----------------------------------------------------------------------------------------------------------
    out=${moduleDir}/${ID}
    nanopore_reads=${clean_readDir}/${ID}_ont.fastq
    bins=${refinedbinDir}/${ID}/metawrap_${min_completion}_${max_contam}_bins
    nanopore=true
    threads=$SLURM_CPUS_PER_TASK
    mem=$Total_Gb

# Checks -------------------------------------------------------------------------------------------------------------------
    if [ ! -s $SOFT/sort_contigs.py ]; then
        error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same folder as the main scripts and all the paths in the config.sh file are correct"
    fi


########################               BEGIN REASSEMBLY PIPELINE!               ########################

module load anaconda3/2023.09
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
conda init
conda activate metawrap-env


STEP="Set up output directory and copy bins"
H1 "$STEP"

    # Make the main output directory
    if [ ! -d $out ]; then
        mkdir $out;
    else
        warning "Warning: $out already exists!"
    fi

    # Reset original_bins folder
    if [ -d ${out}/original_bins ]; then rm -r ${out}/original_bins; fi

    # Copy original bin files
    if [ "$mdmcleaner" = true ]; then
        mkdir ${out}/original_bins
        for i in $bins/*/*_kept_contigs.fasta.gz; do gunzip -k $i && mv ${i%.gz} ${out}/original_bins; done
    else
        cp -r $bins ${out}/original_bins
    fi

    # Make binned_assembly folder
    if [ ! -d ${out}/binned_assembly ]; then mkdir ${out}/binned_assembly; fi

    # Combine the bins into one big assembly file
    if [ ! -s ${out}/binned_assembly/assembly.fa ]; then
        for i in $(ls ${out}/original_bins); do cat ${out}/original_bins/$i >> ${out}/binned_assembly/assembly.fa; done
    fi


########################        RECRUITING READS TO BINS FOR REASSEMBLY         ########################

STEP="Recruiting Reads to Bins"
H1 "$STEP"

    # Set ULIMIT
    ulimit -n 10000
    if [[ $? -ne 0 ]]; then
        ULIMIT=$(ulimit -n)
        warning "Your operating system will allow you to process up to $ULIMIT files at a time. If this is number is less than 4X times of the number of bins you are reassembling, you will likely get an error. Try re-assembling a smaller number of bins."
    fi

    # Make read directory
    mkdir ${out}/reads_for_reassembly

    orig=$(ls ${out}/original_bins/ | wc -l)
    filt=$(ls ${out}/reads_for_reassembly/ | grep .filtered.fastq | wc -l)

    if [ $orig != $filt ]; then
        # Index assembly
        bwa index ${out}/binned_assembly/assembly.fa

        # ONT alignment
        comm "Aligning all reads back to entire assembly and splitting reads into individual fastq files based on their bin membership"
        if [ "$nanopore" = true ]; then
            minimap2 -t $threads -ax map-ont ${out}/binned_assembly/assembly.fa $nanopore_reads \
            | ${SOFT}/filter_nanopore_reads_for_bin_reassembly.py ${out}/original_bins ${out}/reads_for_reassembly
        fi
    fi


########################             REASSEMBLING BINS WITH FLYE              ########################

STEP="Reassembling Bins with Flye"
H1 "$STEP"

    comm "Filtering reads"
    for i in $(ls ${out}/reads_for_reassembly/ | grep .nanopore.fastq); do
        base_name=$(echo ${i} | rev | cut -f 3- -d . | rev)
        if [[ -s ${out}/reads_for_reassembly/${base_name}.filtered.fastq ]]; then
            echo "${base_name} is already filtered."
        else
            python3 ${ScriptPath}/helper_scripts/filter_fastq.py ${out}/reads_for_reassembly/${base_name}.nanopore.fastq ${out}/reads_for_reassembly/${base_name}.filtered.fastq
        fi
    done

    assemble_flye() {
        base_name=$(echo ${1} | rev | cut -f 3- -d . | rev)
        n_reads=${base_name}.filtered.fastq
        bin_name=${base_name}.nanopore
        if [[ -s ${out}/reassemblies/${bin_name}/assembly.fasta ]]; then
            echo -e "\nLooks like $bin_name was already re-assembled. Skipping...\n"
        else
            echo -e "\n--> [ NOW REASSEMBLING ${bin_name} ]\n"

            flye --nano-raw ${out}/reads_for_reassembly/$n_reads --out-dir ${out}/reassemblies/${bin_name} -t $2 --scaffold

            if [[ ! -s ${out}/reassemblies/${bin_name}/assembly.fasta ]]; then
                echo -e "--> XxX Something went wrong with reassembling ${bin_name} XxX"
            else
                echo "--> ${bin_name} was reassembled successfully!"
            fi
        fi
    }

    comm "Reassembly"

    module load flye
    mkdir ${out}/reassemblies

    for i in $(ls ${out}/reads_for_reassembly/ | grep .filtered.fastq); do
        if [ "$nanopore" = true ]; then
            assemble_flye $i $threads
        else
            assemble $i $threads
        fi
    done

    module unload flye
    module load anaconda3/2023.09
    source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
    conda activate metawrap-env

    # Removing short contigs and placing reassemblies in the final folder
    comm "Removing contigs <${len} bp (probably artifacts...)"
    mkdir ${out}/reassembled_bins

    for i in $( ls ${out}/reassemblies/ ); do
        flye_folder=${out}/reassemblies/$i
        bin_name=${flye_folder##*/}

        if [ -s ${out}/reassemblies/${bin_name}/assembly.fasta ]; then

            ${ScriptPath}/helper_scripts/rename_flye_assembly.py \
                ${out}/reassemblies/${bin_name}/assembly.fasta \
                ${out}/reassemblies/${bin_name}/assembly_info.txt \
                > ${out}/reassemblies/${bin_name}/renamed_assembly.fasta

            python3 ${SOFT}/rm_short_contigs.py $len \
             ${out}/reassemblies/${bin_name}/renamed_assembly.fasta \
             > ${out}/reassemblies/${bin_name}/long_scaffolds.fasta

            if [ -s ${out}/reassemblies/${bin_name}/long_scaffolds.fasta ]; then
                echo "$bin_name was reassembled! Processing..."
                mv ${out}/reassemblies/${bin_name}/long_scaffolds.fasta ${out}/reassembled_bins/${bin_name}.fa
            else
                comm "$bin_name was reassembled, but did not yield contigs $len bp. It is possible there were not enough reads."
            fi
        else
            comm "$bin_name was not successfully reassembled. It is possible there were not enough reads."
        fi
    done


########################             RUN CHECKM ON REASSEMBLED BINS             ########################

STEP="CheckM on Reassembled Bins"
H1 "$STEP"

    # Determine --pplacer_threads count. It is either the max thread count or RAM/4, whichever is higher
    ram_max=$(($mem / 40))
    if (( $ram_max < $threads )); then
        p_threads=$ram_max
    else
        p_threads=$threads
    fi
    comm "There is $mem RAM and $threads threads available, and each pplacer thread uses ~40GB, so I will use $p_threads threads for pplacer"

    # Copy over original bins
    for base in $( ls ${out}/original_bins/ | grep "\.fa$" ); do
        i=${out}/original_bins/$base
        cp $i ${out}/reassembled_bins/${base%.*}.orig.fa
    done

    comm "Running CheckM on best bins (reassembled and original)"
    checkm data setRoot $CHECKM_DB
    if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
    mkdir ${out}/tmp
    checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
    if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi

    ${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
    if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
    rm -r ${out}/tmp


########################             FINDING THE BEST VERSION OF EACH BIN             ########################

STEP="Selecting Best Bins"
H1 "$STEP"

    if [ ! -d ${out}/reassembled_best_bins ]; then mkdir ${out}/reassembled_best_bins; fi
    for i in $(${SOFT}/choose_best_bin.py ${out}/reassembled_bins.stats $comp $cont); do
        echo "Copying best bin: $i"
        cp ${out}/reassembled_bins/${i}.fa ${out}/reassembled_best_bins
    done

    o=$(ls -l ${out}/reassembled_best_bins | grep orig | wc -l)
    n=$(ls -l ${out}/reassembled_best_bins | grep nanopore | wc -l)

    announcement "Reassembly results are in! $n bins were improved with reassembly and $o bins were not improved by any reassembly, and thus will stay the same."

    if [[ $(ls ${out}/reassembled_best_bins | wc -l) -gt 0 ]]; then
        comm "Seems that the reassembly went well. You will find the final, best, reassembled bins in ${out}/reassembled_bins, and all intermediate files in ${out}/work_files (which we recommend you delete to save space after you confirm that the pipeline worked)"

        mkdir ${out}/work_files
        mv ${out}/reassembled_bins ${out}/work_files/
        mv ${out}/reassembled_bins.checkm ${out}/work_files/
        # mv ${out}/reassembled_bins.stats ${out}/work_files/
        mv ${out}/reads_for_reassembly ${out}/work_files/
        mv ${out}/nanopore_reads_for_reassembly ${out}/work_files/
        mv ${out}/binned_assembly ${out}/work_files/
        mv ${out}/reassemblies ${out}/work_files/
        rm -r ${out}/original_bins
        mv ${out}/reassembled_best_bins ${out}/reassembled_bins

        Intermediate_files+=(${out}/work_files)
    else
        error "there are no good bins found in ${out}/reassembled_best_bins - something went wrong with choosing the best bins between the reassemblies."
    fi


########################             FINAL CHECKM ON BEST BINS             ########################

STEP="Final CheckM on Best Bins"
H1 "$STEP"

    comm "Re-running CheckM on the best reassembled bins."
    if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
    mkdir ${out}/tmp
    checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
    if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi
    rm -r ${out}/tmp

    comm "Finalizing CheckM stats..."
    ${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
    if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi

    comm "Making CheckM plot of ${out}/reassembled_bins bins"
    checkm bin_qa_plot -x fa ${out}/reassembled_bins.checkm ${out}/reassembled_bins ${out}/reassembled_bins.plot
    if [[ ! -s ${out}/reassembled_bins.plot/bin_qa_plot.png ]]; then warning "Something went wrong with making the CheckM plot. Exiting."; fi
    mv ${out}/reassembled_bins.plot/bin_qa_plot.png ${out}/reassembled_bins.png
    rm -r ${out}/reassembled_bins.plot

    comm "you will find the info on the final reassembled bins in ${out}/reassembled_bins.stats, and a figure summarizing it in ${out}/reassembled_bins.png"

    comm "making reassembly N50, completion, and contamination summary plots."
    head -n 1 ${out}/reassembled_bins.stats > ${out}/original_bins.stats
    grep orig ${out}/reassembled_bins.stats >> ${out}/original_bins.stats
    ${SOFT}/plot_reassembly.py $out $comp $cont ${out}/reassembled_bins.stats ${out}/original_bins.stats
    if [[ $? -ne 0 ]]; then error "Something went wrong with plotting the reassembly summary plots. Exiting..."; fi

    rm -r ${out}/reassembled_bins.checkm

    Intermediate_files+=(${out}/work_files)

    comm "you will find the final bins in ${out}/reassembled_bins"

conda deactivate
module unload anaconda3/2023.09

# Completion
    outputFile="${out}/reassembled_bins.stats"
    step_completion "${outputFile}"

    module_completion

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
