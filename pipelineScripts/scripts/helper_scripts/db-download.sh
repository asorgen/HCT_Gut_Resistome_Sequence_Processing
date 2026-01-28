#!/bin/bash

#SBATCH --partition=DTN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=8GB
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=asorgen@uncc.edu
#SBATCH --job-name db-wget
#SBATCH --time=10:00:00
#SBATCH --output=/scratch/%u/dbLOGs/%x.%j.log

# Set function for output comments -----------------------------------------------------------------------------------------
    H1 () { print_header.py "$1" "H1"; }
    H2 () { print_header.py "$1" "H2"; }
    H3 () { print_header.py "$1" "H3"; }
    comment () { print_header.py "$1" "#"; echo; }
    error () { echo $1; exit 1; }
    pFunc () { echo $1; echo; }


# Set up for script --------------------------------------------------------------------------------------------------------
	H1 "db-download.sh"

	H2 "Job Context"

	comment "Job: $SLURM_JOB_NAME with ID $SLURM_JOB_ID"
	comment "Running on host: `hostname`"

	Total_Gb=$(( SLURM_MEM_PER_NODE / 1024 ))

	JobTime=$(squeue -h -j $SLURM_JOBID -o "%l")
	 
	comment "----- Resources Requested -----"
	comment "Nodes:            $SLURM_NNODES"
	comment "Cores / node:     $SLURM_NTASKS"
	comment "Total memory:     $Total_Gb Gb"
	comment "Wall-clock time:  $JobTime"
	comment "-------------------------------"

# Load environments --------------------------------------------------------------------------------------------------------

	H2 "Modules"
	# module load blast
	# module list


	SECONDS=0

# GTDB-tk database download
	ROOT=/projects/afodor_research3/asorgen/HCT_Gut_Resistome_Pipeline/sequence_processing
	if [[ ! -d ${ROOT}/databases/GTDBtk/release220 ]]; then
	    H2 "GTDB-tk database download"

	    cd ${ROOT}/databases
	    wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
	    tar xvzf gtdbtk_data.tar.gz
	fi

H3 "COMPLETE :)"
duration=$SECONDS
comment "$(elapsed_time "$duration")"
