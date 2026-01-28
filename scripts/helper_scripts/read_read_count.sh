#!/bin/bash

#SBATCH --partition=Orion 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16 
#SBATCH --mem=32GB 
#SBATCH --time=24:00:00 
##SBATCH --job-name=read_count 
##SBATCH --output ${pipeline}/LOGs/read_count.%A.log

# export pipeline=UNC_short
# sbatch --output ${pipeline}/LOGs/read_count.%A.log ./pipelineScripts/scripts/helper_scripts/read_read_count.sh $pipeline 

pipeline=$1
source pipelineScripts/configs/${pipeline}-read.config

echo -e "SampleID\tR1\tR2\tTotal" > ${pipeline}/${pipeline}_tables/${pipeline}_read_counts.tsv

for s in $(tail -n +2 $sampleList); do
    ID="${s%%_*}"
    echo ${ID}

    # Count number of reads per file (1 read = 4 lines in FASTQ)
    num_R1=$(zcat ${pipeline}/0.4_host_decontamination/${ID}_1.fastq.gz | echo -e "$((`wc -l`/4))")
    num_R2=$(zcat ${pipeline}/0.4_host_decontamination/${ID}_2.fastq.gz | echo -e "$((`wc -l`/4))")
    total=$(( num_R1 + num_R2 ))

    echo -e "${ID}\t${num_R1}\t${num_R2}\t${total}" >> ${pipeline}/${pipeline}_tables/${pipeline}_read_counts.tsv

done
