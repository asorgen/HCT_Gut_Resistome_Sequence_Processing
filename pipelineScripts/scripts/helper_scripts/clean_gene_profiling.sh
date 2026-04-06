#!/bin/bash
show_help() {
  cat << EOF
Usage: $(basename "$0") [OPTIONS]

This script removes specific outputs from 5.2_gene_profiling so that the
module can be rerun.

Options:
  -h, --help              Show this help message and exit
  -r, --reset [OPTION]      Optional: Specify the module you would like to reset
                             prodigal
                             bakta
                             amr
                             rgi
                             all
  -p, --pipeline [OPTION]   Specify the dataset you would like to reset
                             Duke_short
                             Duke_long
                             Duke_hybrid
                             UNC_short
  -s, --sequences [OPTION]  Specify the sequence type to reset
                             asm
                             bin
Example:
  $(basename "$0") -r rgi -p UNC_short -s asm
EOF
}

# Set up -----------------------------------------------------------------------------------------------------------
	# Initialize variables
	pipelines=()
	resets=()
	sequences=()
	r_option=""
	p_option=""
	verbose=false

	# Parse options
		while [[ "$#" -gt 0 ]]; do
		  case "$1" in
		    -h|--help)
		      show_help
		      exit 0
		      ;;
		    -r|--reset)
		    	shift
		    	while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
		    		resets+=("$1")
		    		shift
		    	done
		      ;;
		    -p|--pipeline)
		    	shift
		    	while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
		    		pipelines+=("$1")
		    		shift
		    	done
		      ;;
		    -s|--sequences)
		    	shift
		    	while [[ "$#" -gt 0 && ! "$1" =~ ^- ]]; do
		    		sequences+=("$1")
		    		shift
		    	done
		      ;;
		    -v|--verbose)
		      verbose=true
		      shift
		      ;;
		    -*)
		      echo "Unknown option: $1"
		      show_help
		      exit 1
		      ;;
		    *)
		      echo "Unexpected argument: $1"
		      show_help
		      exit 1
		      ;;
		  esac
		done

	# Validate required input
	if [ "${#pipelines[@]}" -eq 0 ]; then
	  echo "Error: At least one pipeline must be specified with -p or --pipeline."
	  show_help
	  exit 1
	fi


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../configs/private.config"
root=${HPC_PROJECTS}/HCT_Gut_Resistome_Study/HCT_Gut_Resistome_Pipeline/sequence_processing
# pipelines=(Duke_short Duke_long Duke_hybrid UNC_short)

for pipeline in "${pipelines[@]}"; do
	echo $pipeline
	
	for clean in "${resets[@]}"; do 
		echo $clean
		
		for method in "${sequences[@]}"; do
			echo $method
			if [[ $method == 'asm' ]]; then
				module_directory=5.5_AA_amr_assembly
			else
				module_directory=5.6_AA_amr_bins
			fi

			dir=${root}/${pipeline}/${module_directory}
			if [[ ! -d $dir ]]; then continue; fi
			echo $dir

			sampleIDs=$(ls ${dir})

			for ID in ${sampleIDs}; do
				
				if [[ $ID == "logs" || $ID == "COMPLETE" ]]; then continue; fi
			if [[ -n "${exclude_ids}" ]] && echo "${exclude_ids}" | grep -qw "$ID"; then continue; fi
				
				prodigalFile=${dir}/${ID}/${ID}_genes.faa
				baktaDir=${dir}/${ID}/bakta
				amrFile=${dir}/${ID}/${ID}.amrfinder.tsv
				rgiFiles=$(ls ${dir}/${ID}/*rgi*)
				# echo "${rgiFiles[@]}"
				genelistFile=${dir}/${ID}/${ID}_geneslist.tsv
				annotationFile=${dir}/${ID}/${ID}_gene_annotations.tsv

				if [[ $clean == "prodigal" && -f $prodigalFile ]]; then
					rm $prodigalFile
					rm $genelistFile
				fi

				if [[ $clean == "bakta" && -d $baktaDir ]]; then
					rm -r $baktaDir
				fi

				if [[ $clean == "amr" && -f $amrFile ]]; then
					rm $amrFile
				fi

				if [[ $clean == "rgi" && -n "$rgiFiles" ]]; then
					# echo "Files found"
					rm $rgiFiles
				fi

				if [[ $clean == "all" ]]; then
					rm -r ${dir}/${ID}
				fi

				if [[ -f $annotationFile ]]; then
					rm $annotationFile
				fi
			done # for ID in ${sampleIDs}
		done # for method in asm bin
	done # for clean in "${resets[@]}"
done # for pipeline in "${pipelines[@]}"


