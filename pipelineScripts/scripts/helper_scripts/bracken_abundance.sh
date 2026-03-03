#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=24gb
#SBATCH --time=02:00:00
#SBATCH --job-name=brackenJoin
#SBATCH --output=LOGS/bracken-join.%A.log


brackenDir=$1
readType=$2

python3 ${HOME}/PROGRAMS/Bracken-2.7/analysis_scripts/combine_bracken_outputs.py \
--files ${brackenDir}/*.out \
-o ${brackenDir}/${readType}_bracken_counts.tsv
