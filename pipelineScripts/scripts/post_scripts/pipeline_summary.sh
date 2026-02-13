#!/bin/bash

# sbatch --output=${dataset}_summary.%A.log pipeline_summary.sh -p ${dataset}

#SBATCH --partition=Orion
#SBATCH --job-name=pipeline_summary
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=8GB
#SBATCH --time=01:00:00

source /users/asorgen/.bashrc
config_file=$(which config-metawrap)
source $config_file
# source short-read.config

# --- Defaults ---
PIPELINES=()
MODULES=()

# --- Help message ---
show_help() {
  cat << EOF

Usage: $(basename "$0") [-p PIPELINES] [-m MODULES]

This script summarizes all output files from pipeline modules across samples into single, individual files.

Options:
  -h, --help       Show this help message and exit

  Required:
    -p, --pipeline Comma-separated list of pipeline names.
                    Valid options:
                     - Duke_short
                     - Duke_long
                     - Duke_hybrid
                     - UNC_short

  Optional:
    -m, --modules  Comma-separated list of module names.
                    Valid options:
                     - all (default)
                     - read_qc
                     - evaluation
                     - kraken2
                     - metaphlan4
                     - binning
                     - refine_bins
                     - reassemble_bins
                     - classify_bins
                     - annotate_bins
                     - amr_detection
                     - asm_gene_profiling
                     - bin_gene_profiling
                     - shortbred
                     - rgi_bwt

Examples:
  $(basename "$0") -p UNC_short -m kraken2,shortbred
  sbatch $(basename "$0") -p Duke_short 
EOF
}

# --- Parse arguments ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -p|--pipeline)
      IFS=',' read -r -a PIPELINES <<< "$2"
      shift 2
      ;;
    -m|--module)
      IFS=',' read -r -a MODULES <<< "$2"
      shift 2
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      show_help
      exit 1
      ;;
  esac
done

module load anaconda3/2023.09
source /apps/pkg/anaconda3/2023.09/etc/profile.d/conda.sh
conda activate metawrap-env

# Set function for output comments
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo "> $1"; echo; }

if [[ -v SLURM_JOB_ID ]]; then
	H1 "Job Context"
	    export OMP_NUM_THREADS=$SLURM_NTASKS
	    comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
	    comment "Running on host: `hostname`"

	    Total_Gb=$(( SLURM_MEM_PER_NODE / 1000 ))

	    JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")

	    echo 
	    comment "----- Resources Requested -----"
	    comment "Nodes:            $SLURM_NNODES"
	    comment "Cores / node:     $SLURM_NTASKS"
	    comment "Total memory:     $Total_Gb Gb"
	    comment "Wall-clock time:  $JobTime"
	    comment "-------------------------------"
fi


# Module defaults
read_qc=false
evaluation=false
kraken2=false; metaphlan4=false
binning=false; refine_bins=false; reassemble_bins=false
classify_bins=false; annotate_bins=false
amr_detection=false; asm_gene_profiling=false; bin_gene_profiling=false; shortbred=false; rgi_bwt=false
all=true

if [[ ${#MODULES[@]} -gt 0 ]]; then
	all=false
	for module in "${MODULES[@]}"; do
	    if [[ "$module" == "read_qc" || "$module" == "all" ]]; then read_qc=true; fi
	    if [[ "$module" == "evaluation" || "$module" == "all" ]]; then evaluation=true; fi
	    if [[ "$module" == "kraken2" || "$module" == "all" ]]; then kraken2=true; fi
	    if [[ "$module" == "metaphlan4" || "$module" == "all" ]]; then metaphlan4=true; fi
	    if [[ "$module" == "binning" || "$module" == "all" ]]; then binning=true; fi
	    if [[ "$module" == "refine_bins" || "$module" == "all" ]]; then refine_bins=true; fi
	    if [[ "$module" == "reassemble_bins" || "$module" == "all" ]]; then reassemble_bins=true; fi
	    if [[ "$module" == "classify_bins" || "$module" == "all" ]]; then classify_bins=true; fi
	    if [[ "$module" == "annotate_bins" || "$module" == "all" ]]; then annotate_bins=true; fi
	    if [[ "$module" == "amr_detection" || "$module" == "all" ]]; then amr_detection=true; fi
	    if [[ "$module" == "shortbred" || "$module" == "all" ]]; then shortbred=true; fi
	    if [[ "$module" == "rgi_bwt" || "$module" == "all" ]]; then rgi_bwt=true; fi
	    if [[ "$module" == "amr_rgi_merge" || "$module" == "all" ]]; then amr_rgi_merge=true; fi
	    if [[ "$module" == "asm_gene_profiling" || "$module" == "all" ]]; then asm_gene_profiling=true; fi
	    if [[ "$module" == "bin_gene_profiling" || "$module" == "all" ]]; then bin_gene_profiling=true; fi
	    if [[ "$module" == "all" ]]; then all=true; fi
	done
fi

# If all modules are requested, enable every module flag
if $all; then
	read_qc=true; evaluation=true
	kraken2=true; metaphlan4=true
	binning=true; refine_bins=true; reassemble_bins=true
	classify_bins=true; annotate_bins=true
	amr_detection=true; asm_gene_profiling=true; bin_gene_profiling=true; shortbred=true; rgi_bwt=true
fi

for pipeline in "${PIPELINES[@]}"; do
	H1 "$pipeline"

	source pipelineScripts/configs/${pipeline}-read.config
	export ps_path=${ROOT}/pipelineScripts/scripts/post_scripts
	# echo "post scripts: $ps_path"

	cd $datasetDir

	out=${PROCESSED_ROOT}/${pipeline}_tables
	mkdir -p ${out}

	# read_qc
		
		if $read_qc; then 
			if [[ $pipeline = "Duke_short" || $pipeline = "Duke_long" || $pipeline = "UNC_short" ]]; then
				SECONDS=0
				H2 "read_qc"
				module load multiqc

				module_dir=0.1_pre_qc
					if [[ -d ${module_dir}/multiqc_data ]]; then rm -r ${module_dir}/multiqc_data; rm ${module_dir}/multiqc_report.html; fi
					multiqc ${module_dir}/ -o ${module_dir}/
					if [[ $? -ne 0 ]]; then error "Something went wrong with the pre-qc report. Exiting..."; fi
					mv ${module_dir}/multiqc_data/multiqc_general_stats.txt ${out}/${pipeline}_pre-QC_report.tsv

				module_dir=0.2_deduplication
					if [[ -d ${module_dir}/multiqc_data ]]; then rm -r ${module_dir}/multiqc_data; rm ${module_dir}/multiqc_report.html; fi
					multiqc ${module_dir}/QC_report -o ${module_dir}/
					if [[ $? -ne 0 ]]; then error "Something went wrong with the pre-qc report. Exiting..."; fi
					mv ${module_dir}/multiqc_data/multiqc_general_stats.txt ${out}/${pipeline}_dedup_report.tsv

				module_dir=0.3_sequence_trim
					if [[ -d ${module_dir}/multiqc_data ]]; then rm -r ${module_dir}/multiqc_data; rm ${module_dir}/multiqc_report.html; fi
					multiqc ${module_dir}/QC_report -o ${module_dir}/
					if [[ $? -ne 0 ]]; then error "Something went wrong with the pre-qc report. Exiting..."; fi
					mv ${module_dir}/multiqc_data/multiqc_general_stats.txt ${out}/${pipeline}_dedup_trim_report.tsv

				module_dir=0.4_host_decontamination
					if [[ -d ${module_dir}/multiqc_data ]]; then rm -r ${module_dir}/multiqc_data; rm ${module_dir}/multiqc_report.html; fi
					multiqc ${module_dir}/QC_report -o ${module_dir}/
					if [[ $? -ne 0 ]]; then error "Something went wrong with the pre-qc report. Exiting..."; fi
					mv ${module_dir}/multiqc_data/multiqc_general_stats.txt ${out}/${pipeline}_post-QC_report.tsv

				module unload multiqc
				comment "COMPLETE :)"
				duration=$SECONDS
				comment "$(elapsed_time "$duration")"
				# if ! $all; then exit 0; fi
			fi
		fi

	# evaluation
		module_dir=1.2_evaluation
		if $evaluation; then 
			SECONDS=0
			H2 "evaluation"
			
		  # if [[ -s ${out}/${pipeline}_draft_assembly_stats.tsv ]]; then rm ${out}/${pipeline}_draft_assembly_stats.tsv; fi
			# statswrapper.sh ${module_dir}/*_draft_assembly.fasta >> ${out}/${pipeline}_draft_assembly_stats.tsv
			# if [[ $? -ne 0 ]]; then error "Something went wrong with the draft assembly statswrapper. Exiting..."; fi

		  # if [[ -s ${out}/${pipeline}_final_assembly_stats.tsv ]]; then rm ${out}/${pipeline}_final_assembly_stats.tsv; fi
			# statswrapper.sh ${module_dir}/*_final_assembly.fasta >> ${out}/${pipeline}_final_assembly_stats.tsv
			# if [[ $? -ne 0 ]]; then error "Something went wrong with the final assembly statswrapper. Exiting..."; fi

		  # if [[ -d ${module_dir}/QUAST_out ]]; then rm -r ${module_dir}/QUAST_out; fi
		  # mkdir -p ${module_dir}/QUAST_out
		  # module load quast
		  # quast.py -t $SLURM_NTASKS -o ${module_dir}/QUAST_out -m 500 ${module_dir}/*.fasta
		  # if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi
		  # module unload quast

		  # cd ${module_dir}/QUAST_out/
		  # for file in *; do mv "$file" "${pipeline}_$file"; done        
		  # cd -
		  # mv ${module_dir}/QUAST_out/${pipeline}_report.tsv ${out}/${pipeline}_QUAST_report.tsv
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# kraken2
		module_dir=2.1_kraken2
		module_dir=2.2_bracken
		if $kraken2; then 
			SECONDS=0
			H2 "kraken2"

		  # Assembly output
		  if [[ -d ${module_dir}/assembly ]]; then
			  H3 "Assemblies"
			  if [[ -s ${out}/${pipeline}_assembly_bracken_counts.tsv ]]; then rm ${out}/${pipeline}_assembly_bracken_counts.tsv; fi
			  python3 ${ps_path}/combine_bracken_outputs.py --files ${module_dir}/assembly/*.out -o ${out}/${pipeline}_assembly_bracken_counts.tsv
			  if [[ $? -ne 0 ]]; then error "Something went wrong while compiling the Bracken count table. Exiting..."; fi
		  fi

			# Short-read ouputs
			if [[ -d ${module_dir}/sr ]]; then
				H3 "Short reads"
				if [[ -s ${out}/${pipeline}_sr_bracken_counts.tsv ]]; then rm ${out}/${pipeline}_sr_bracken_counts.tsv; fi
				python3 ${ps_path}/combine_bracken_outputs.py --files ${module_dir}/sr/*.out -o ${out}/${pipeline}_sr_bracken_counts.tsv
				if [[ $? -ne 0 ]]; then error "Something went wrong while compiling the Bracken count table. Exiting..."; fi

			fi

			# Long-read outputs
			if [[ -d ${module_dir}/ont ]]; then
					H3 "Long reads"
			    if [[ -s ${out}/${pipeline}_ont_bracken_counts.tsv ]]; then rm ${out}/${pipeline}_ont_bracken_counts.tsv; fi
			    python3 ${ps_path}/combine_bracken_outputs.py --files ${module_dir}/ont/*.out -o ${out}/${pipeline}_ont_bracken_counts.tsv
			    if [[ $? -ne 0 ]]; then error "Something went wrong while compiling the Bracken count table. Exiting..."; fi

			fi

			comment "COMPLETE :)"
			duration=$SECONDS
			comment "$(elapsed_time "$duration")"
			# if ! $all; then exit 0; fi
		fi

	# binning
		module_dir=3.1_binning
		if $binning; then 
			SECONDS=0
			H2 "binning"

			# echo -e "Pipeline\tSampleID\tBinner\tBins" > ${out}/${pipeline}_initial_bin_summary.tsv
			# for dir in `ls ${module_dir}/`; do
			  
			#   if [[ -d ${module_dir}/${dir}/metabat2_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/metabat2_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tmetabat2\t${count}" >> ${out}/${pipeline}_initial_bin_summary.tsv
			#   fi

			#   if [[ -d ${module_dir}/${dir}/maxbin2_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/maxbin2_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tmaxbin2\t${count}" >> ${out}/${pipeline}_initial_bin_summary.tsv
			#   fi

			#   if [[ -d ${module_dir}/${dir}/concoct_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/concoct_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tconcoct\t${count}" >> ${out}/${pipeline}_initial_bin_summary.tsv
			#   fi

			# done
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# refine_bins
		module_dir=3.2_refine_bins
		if $refine_bins; then 
			SECONDS=0
			H2 "refine_bins"

			# echo -e "Pipeline\tSampleID\tBinner\tBins" > ${out}/${pipeline}_refined_bin_summary.tsv
			# for dir in `ls ${module_dir}/`; do
			  
			#   if [[ -d ${module_dir}/${dir}/metabat2_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/metabat2_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tmetabat2\t${count}" >> ${out}/${pipeline}_refined_bin_summary.tsv
			#   fi

			#   if [[ -d ${module_dir}/${dir}/maxbin2_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/maxbin2_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tmaxbin2\t${count}" >> ${out}/${pipeline}_refined_bin_summary.tsv
			#   fi

			#   if [[ -d ${module_dir}/${dir}/concoct_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/concoct_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tconcoct\t${count}" >> ${out}/${pipeline}_refined_bin_summary.tsv
			#   fi

			# if [[ -d ${module_dir}/${dir}/metawrap_70_10_bins ]]; then
			#   	count=$(ls ${module_dir}/${dir}/metawrap_70_10_bins | wc -l)
			#   	echo -e "${pipeline}\t${dir}\tmetawrap\t${count}" >> ${out}/${pipeline}_refined_bin_summary.tsv
			#   fi

			# done

			# if [[ -f ${out}/${pipeline}_refined_bin_stats.tsv ]]; then rm ${out}/${pipeline}_refined_bin_stats.tsv; fi
			# ${ps_path}/refine_stats.sh ${module_dir} metabat2_bins.stats ${pipeline}_refined_bin_stats.tsv	
			# ${ps_path}/refine_stats.sh ${module_dir} maxbin2_bins.stats ${pipeline}_refined_bin_stats.tsv
			# ${ps_path}/refine_stats.sh ${module_dir} concoct_bins.stats ${pipeline}_refined_bin_stats.tsv
			# ${ps_path}/refine_stats.sh ${module_dir} metawrap_70_10_bins.stats ${pipeline}_refined_bin_stats.tsv

			# mv ${module_dir}/${pipeline}_refined_bin_stats.tsv ${out}/${pipeline}_refined_bin_stats.tsv
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# reassemble_bins
		module_dir=3.3_reassemble_bins
		if $reassemble_bins; then 
			SECONDS=0
			H2 "reassemble_bins"

			# if [[ -f ${out}/${pipeline}_reassembled_bin_stats.tsv ]]; then rm ${out}/${pipeline}_reassembled_bin_stats.tsv; fi
			# # ${ps_path}/reassem_stats.sh reassemble_bins original_bins.stats ${pipeline}_reassembled_bin_stats.tsv
			# ${ps_path}/reassem_stats.sh ${module_dir} reassembled_bins.stats ${pipeline}_reassembled_bin_stats.tsv
			# mv ${module_dir}/${pipeline}_reassembled_bin_stats.tsv ${out}/${pipeline}_reassembled_bin_stats.tsv
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# classify_bins
		module_dir=4.1_classify_bins
		if $classify_bins; then 
			SECONDS=0
			H2 "classify_bins"

			# if [[ -f ${out}/${pipeline}_bin_classification.tsv ]]; then rm ${out}/${pipeline}_bin_classification.tsv; fi
			# ${ps_path}/classify_bins_summary.sh ${module_dir} bin_taxonomy.tab ${pipeline}_bin_classification.tsv
			# mv ${module_dir}/${pipeline}_bin_classification.tsv ${out}/${pipeline}_bin_classification.tsv
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# annotate_bins
		module_dir=4.2_annotate_bins
		if $annotate_bins; then 
			SECONDS=0
			H2 "annotate_bins"

			# if [[ -f ${out}/${pipeline}_bin_annotation.tsv ]]; then rm ${out}/${pipeline}_bin_annotation.tsv; fi
			# ${ps_path}/annotate_bins_summary.sh ${module_dir} bin_funct_annotations ${pipeline}_bin_annotation.tsv
			# mv ${module_dir}/${pipeline}_bin_annotation.tsv ${out}/${pipeline}_bin_annotation.tsv
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# amr_detection
		module_dir=5.1_amr_detection
		if $amr_detection; then 
			SECONDS=0
		    H2 "amr_detection"

		  #   H2 "RGI"
		  #   cd ${module_dir}/RGI
		  #   python3 ${ps_path}/rgi_count.py -f . -s $pipeline
		  #   if [[ $? -ne 0 ]]; then error "Something went wrong with RGI. Exiting..."; fi
		  #   cd -
			# 	mv ${module_dir}/RGI/RGI_${pipeline}.tsv ${out}/${pipeline}_RGI_counts.tsv

		  #   H2 "AMRFinder"
		  #   cd ${module_dir}/AMR
		  #   python3 ${ps_path}/amr_count.py -f . -s $pipeline
		  #   if [[ $? -ne 0 ]]; then error "Something went wrong with AMRFinder. Exiting..."; fi
		  #   cd -
			# mv ${module_dir}/AMR/AMRFinder_${pipeline}.tsv ${out}/${pipeline}_AMRFinder_counts.tsv
			# comment "COMPLETE :)"
			# duration=$SECONDS
			# comment "$(elapsed_time "$duration")"
			# # if ! $all; then exit 0; fi
		fi

	# asm_gene_profiling
		module_dir=5.2_gene_profiling/asm
		if $asm_gene_profiling; then 
		  SECONDS=0
		  H2 "Assembly Gene Profiling"

		#   if [[ -f ${out}/${pipeline}_asm_gene_annotations.tsv ]]; then rm ${out}/${pipeline}_asm_gene_annotations.tsv; fi
		  
		#   ${ps_path}/gene_profiling_summary.sh \
		# 	  ${module_dir} \
		# 	  gene_annotations.tsv \
		# 	  ${pipeline}_asm_gene_annotations.tsv
			  
		#   mv ${module_dir}/${pipeline}_asm_gene_annotations.tsv ${out}/${pipeline}_asm_gene_annotations.tsv
		#   comment "COMPLETE :)"
		#   duration=$SECONDS
		#   comment "$(elapsed_time "$duration")"
		#   if ! $all; then exit 0; fi
		# fi

		# module_dir=5.2_gene_profiling/bin
		# if $bin_gene_profiling; then 
		#   SECONDS=0
		#   H1 "Bin Gene Profiling"
		#   if [[ -f ${out}/${pipeline}_bin_gene_annotations.tsv ]]; then rm ${out}/${pipeline}_bin_gene_annotations.tsv; fi
		  
		#   ${ps_path}/gene_profiling_summary.sh \
		# 	  ${module_dir} \
		# 	  gene_annotations.tsv \
		# 	  ${pipeline}_bin_gene_annotations.tsv
			  
		#   mv ${module_dir}/${pipeline}_bin_gene_annotations.tsv ${out}/${pipeline}_bin_gene_annotations.tsv
		#   comment "COMPLETE :)"
		#   duration=$SECONDS
		#   comment "$(elapsed_time "$duration")"
		#   # if ! $all; then exit 0; fi
		fi

	# shortbred
		module_dir=5.3_shortbred
		if $shortbred; then 
			SECONDS=0
		  H2 "shortbred"

		  python3 ${ps_path}/summarize_shortbred.py -f $module_dir -l $pipeline -r ${out}/${pipeline}_post-QC_report.tsv
		  if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi

			mv ${module_dir}/${pipeline}_shortbred.tsv ${out}/

			comment "COMPLETE :)"
			duration=$SECONDS
			comment "$(elapsed_time "$duration")"
			# if ! $all; then exit 0; fi
		fi

	# rgi_bwt
		module_dir=5.4_rgi_bwt
		if $rgi_bwt; then 
			SECONDS=0
		  H2 "RGI bwt"

		  python3 ${ps_path}/summarize_rgi_bwt.py \
		  -f $module_dir \
		  -l $pipeline \
		  --mapped 5 --mapq 0 --coverage 0 \
		  -r ${out}/${pipeline}_post-QC_report.tsv

		  if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi

			mv ${module_dir}/${pipeline}_rgi_bwt_RPKM.tsv ${out}/

			comment "COMPLETE :)"
			duration=$SECONDS
			comment "$(elapsed_time "$duration")"
			# if ! $all; then exit 0; fi
		fi

	# metaphlan4
		module_dir=2.3_metaphlan4
		if $metaphlan4; then
			SECONDS=0
			H2 "MetaPhlAn4"

		  ${ps_path}/metaphlanMerge.sh -i $module_dir -o ${out}/${pipeline}
		  if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi

		  ${ps_path}/BuildTaxaTable_MetaPhlAn4.sh -i ${out}/${pipeline}_metaphlan4_rel_abun.tsv -o ${out}/${pipeline} -p $module_dir
		  if [[ $? -ne 0 ]]; then error "Something went wrong. Exiting..."; fi


			comment "COMPLETE :)"
			duration=$SECONDS
			comment "$(elapsed_time "$duration")"
			# if ! $all; then exit 0; fi
		fi

	cd -
done


