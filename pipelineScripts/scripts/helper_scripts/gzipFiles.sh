#!/bin/bash

#SBATCH --partition=Orion
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --time=50:00:00
##SBATCH --job-name=gzip
##SBATCH --output=LOGS/gzipFiles.%A.log

# Usage:
# sbatch gzipFiles.sh HCT_Gut_Resistome_Pipeline/sequence_processing/SHORT/clean_reads fastq
# sbatch --job-name=UNC-gzip --output=HCT_Gut_Resistome_Pipeline/sequence_processing/UNC-BMT/LOGs/gzipFiles.%A.log HCT_Gut_Resistome_Pipeline/sequence_processing/SHORT/scripts/helper_scripts/gzipFiles.sh $uncFolder fastq
# sbatch --job-name=short-gzip --output=HCT_Gut_Resistome_Pipeline/sequence_processing/SHORT/LOGs/gzipFiles.%A.log HCT_Gut_Resistome_Pipeline/sequence_processing/SHORT/scripts/helper_scripts/gzipFiles.sh $folder fastq

directory=$1
fileType=$2

# gzip ${directory}/*.${fileType}
count=1
for file in `ls ${directory}`; do
	ID="${file%%.*}"
	
	if [[ -s ${directory}/${ID}.${fileType} ]]; then
				
		if [[ -s ${directory}/${ID}.${fileType}.gz ]]; then
			rm ${directory}/${ID}.${fileType}.gz
		fi
		
		gzip ${directory}/${file}
		echo "${count}. ${file}.gz"
	else
		echo "${count}. ${file}"	
	fi
	count=$((count + 1))

done