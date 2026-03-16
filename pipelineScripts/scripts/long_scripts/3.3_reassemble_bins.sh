#!/bin/bash

##SBATCH --partition=Orion
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=48
##SBATCH --mem=1200GB
##SBATCH --time=48:00:00
##SBATCH --mail-user=${email}
##SBATCH --mail-type=END,FAIL
##SBATCH --job-name=test
##SBATCH --output=${moduleDir}/${ID}/testing.log

# H1 "Job Context"
    OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
    # comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
    # comment "Running on host: `hostname`"

    Total_Gb=$(( SLURM_MEM_PER_NODE / 1000 ))

    JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

    echo 
    echo "----- Resources Requested -----"
    echo "Nodes:            $SLURM_NNODES"
    echo "Cores / node:     $SLURM_CPUS_PER_TASK"
    echo "Total memory:     $Total_Gb Gb"
    echo "Wall-clock time:  $JobTime"
    echo "-------------------------------"

# Set function for output comments
	H1 () { print_header.py "$1" "H1"; }
	H2 () { print_header.py "$1" "H2"; }
	H3 () { print_header.py "$1" "H3"; }
	comment () { print_header.py "$1" "#"; }


# cd /projects/afodor_research3/asorgen/HCT_Gut_Resistome_Pipeline/sequence_processing/LONG

module load anaconda3/2023.09
source ~/.bashrc
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
conda init
conda activate metawrap-env


comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }

# these functions are for parallelizing the reassembly
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

########################               LOADING IN THE PARAMETERS                ########################


	# setting scripts and databases from config file (should be in same folder as main script)
		source $pipelineConfig
		source $config_file
		source $bashrc
		source $bash_profile


	# default params
		threads=1; mem=40; comp=70; cont=10; len=500
		bins=None; f_reads=None; r_reads=None; out=None

	# long options defaults
		strict_max=2; permissive_max=5
		run_checkm=true
		run_parallel=false
		nanopore=false
		mdmcleaner=false

	# load in params
		out=${moduleDir}/${ID}
		nanopore_reads=${clean_readDir}/${ID}_ont.fastq
		bins=${refinedbinDir}/${ID}/metawrap_${max_completion}_${min_contam}_bins
		nanopore=true
		threads=$SLURM_CPUS_PER_TASK
		mem=$Total_Gb

########################           MAKING SURE EVERYTHING IS SET UP             ########################

	# Checks for correctly configures meta-scripts folder
		if [ ! -s $SOFT/sort_contigs.py ]; then
			error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
		fi


########################               BEGIN REASSEMBLY PIPELINE!               ########################

	SECONDS=0
	
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
	announcement "RECRUITING READS TO BINS FOR REASSEMBLY"

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

	announcement "REASSEMBLING BINS WITH FLYE"
	mkdir ${out}/reassemblies
	
	comm "Filtering reads"
	for i in $(ls ${out}/reads_for_reassembly/ | grep .nanopore.fastq); do
		base_name=$(echo ${i} | rev | cut -f 3- -d . | rev)
		if [[ -s ${out}/reads_for_reassembly/${base_name}.filtered.fastq ]]; then
			echo "${base_name} is already filtered."
		else
			python scripts/helper_scripts/filter_fastq.py ${out}/reads_for_reassembly/${base_name}.nanopore.fastq ${out}/reads_for_reassembly/${base_name}.filtered.fastq
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

	for i in $(ls ${out}/reads_for_reassembly/ | grep .filtered.fastq); do

		if [ "$nanopore" = true ]; then
			assemble_flye $i $threads 
		else
			assemble $i $threads 
		fi
	done

	module unload flye

	# removing short contigs and placing reassemblies in the final folder
	conda activate metawrap-env

	comm "Removing contigs <${len} bp (probably artifacts...)"
	mkdir ${out}/reassembled_bins

	for i in $( ls ${out}/reassemblies/ ); do
		flye_folder=${out}/reassemblies/$i #${moduleDir}/D21309D29/reassemblies/bin.9.nanopore
		bin_name=${flye_folder##*/} # bin.9.nanopore		
		
		#remove shortest contigs (probably artifacts...)
		if [ -s ${out}/reassemblies/${bin_name}/assembly.fasta ]; then

			scripts/helper_scripts/rename_flye_assembly.py \
				${out}/reassemblies/${bin_name}/assembly.fasta \
				${out}/reassemblies/${bin_name}/assembly_info.txt \
				> ${out}/reassemblies/${bin_name}/renamed_assembly.fasta

			${SOFT}/rm_short_contigs.py $len \
			 ${out}/reassemblies/${bin_name}/renamed_assembly.fasta\
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
	# determine --pplacer_threads count. It is either the max thread count or RAM/4, whichever is higher
	ram_max=$(($mem / 40))
	if (( $ram_max < $threads )); then
		p_threads=$ram_max
	else
		p_threads=$threads
	fi
	comm "There is $mem RAM and $threads threads available, and each pplacer thread uses ~40GB, so I will use $p_threads threads for pplacer"

	# copy over original bins
	for base in $( ls ${out}/original_bins/ | grep "\.fa$" ); do 
		i=${out}/original_bins/$base
		cp $i ${out}/reassembled_bins/${base%.*}.orig.fa
	done


	comm "Running CheckM on best bins (reassembled and original)"
	if [[ -d ${out}/reassembled_bins.checkm ]]; then rm -r ${out}/reassembled_bins.checkm; fi
	mkdir ${out}/tmp
	checkm lineage_wf -x fa ${out}/reassembled_bins ${out}/reassembled_bins.checkm -t $threads --tmpdir ${out}/tmp --pplacer_threads $p_threads
	if [[ ! -s ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv ]]; then error "Something went wrong with running CheckM. Exiting..."; fi

	${SOFT}/summarize_checkm.py ${out}/reassembled_bins.checkm/storage/bin_stats_ext.tsv | (read -r; printf "%s\n" "$REPLY"; sort) > ${out}/reassembled_bins.stats
	if [[ $? -ne 0 ]]; then error "Cannot make checkm summary file. Exiting."; fi
	rm -r ${out}/tmp


	announcement "FINDING THE BEST VERSION OF EACH BIN"
	if [ ! -d ${out}/reassembled_best_bins ]; then mkdir ${out}/reassembled_best_bins; fi
	for i in $(${SOFT}/choose_best_bin.py ${out}/reassembled_bins.stats $comp $cont); do 
		echo "Copying best bin: $i"
		cp ${out}/reassembled_bins/${i}.fa ${out}/reassembled_best_bins 
	done

	o=$(ls -l ${out}/reassembled_best_bins | grep orig | wc -l)
	n=$(ls -l ${out}/reassembled_best_bins | grep nanopore | wc -l)

	announcement "Reassembly results are in! $n bins were improved with reassembly and $o bins were not improved by any reassembly, and thus will stay the same."

	if [[ $(ls ${out}/reassembled_best_bins | wc -l) -gt 0 ]]; then 
		comm "Seems that the reassembly went well. You will find the final, best, reassembled bins in ${out}/reassembled_bins, and all intermediate files in ${out}/work_files (which we recomend you delete to save space after you confirm that the pipeline worked)"
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
	else
		error "there are no good bins found in ${out}/reassembled_best_bins - something went wrong with choosing the best bins between the reassemblies."
	fi

	comm "Re-running CheckM on the best reasembled bins."
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

	comm "making reassembly N50, compleiton, and contamination summary plots."
	head -n 1 ${out}/reassembled_bins.stats > ${out}/original_bins.stats
	grep orig ${out}/reassembled_bins.stats >> ${out}/original_bins.stats
	${SOFT}/plot_reassembly.py $out $comp $cont ${out}/reassembled_bins.stats ${out}/original_bins.stats
	if [[ $? -ne 0 ]]; then error "Something went wrong with plotting the reassembly summary plots. Exiting..."; fi

rm -r ${out}/reassembled_bins.checkm
rm -r ${out}/work_files

comm "you will find the final bins in ${out}/reassembled_bins"

conda deactivate
module unload anaconda3/2023.09

# Completion status
    if [[ -s ${moduleDir}/${ID}/reassembled_bins.stats ]]; then
        touch ${moduleDir}/${ID}/COMPLETE
    fi

H1 "PIPELINE COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
