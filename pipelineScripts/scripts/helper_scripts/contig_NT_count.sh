#!/bin/bash

#SBATCH --partition=Orion 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=32GB 
#SBATCH --time=24:00:00 
##SBATCH --job-name=read_count 
##SBATCH --output ${pipeline}/LOGs/read_count.%A.log

# pipeline=UNC_short
# sbatch --job-name=${pipeline}_NT_count --output ${pipeline}/LOGs/NT_count.%A.log ./pipelineScripts/scripts/helper_scripts/contig_NT_count.sh $pipeline 

pipeline=$1
source pipelineScripts/configs/${pipeline}-read.config

echo -e "SampleID\tNucleotides" > ${pipeline}/1.2_evaluation/${pipeline}_contig_NT_counts.tsv

for s in $(tail -n +2 $sampleList); do
    ID="${s%%_*}"
    echo ${ID}

    # Count number of nucleotides
    stats.sh in="${pipeline}/1.2_evaluation/${ID}_final_assembly.fasta" \
    | awk -v id="$ID" -F'\t' 'BEGIN{IGNORECASE=1} /^ *All[[:space:]]/ {
        gsub(/,/, "", $5);
        print id "\t" $5; 
        exit
    }' >> "${pipeline}/1.2_evaluation/${pipeline}_contig_NT_counts.tsv"
    
done
