# Taxonomic classification with Kraken2

**Description**: Running Kraken2 on the reads will give us an idea of the taxonomic composition of the community. Running Kraken2 on the assembly will give us an idea what taxonomic groups were assembled better than others (the assembly process is heavily biased and should no be used to infer overall community composition).

## Module overview

  1. Runs Kraken2 on reads
  2. Runs Kraken2 on the assembly
  3. Runs Kraken-translate on the output
  4. Generates a kronagram of all files

## Setup

See [setup](manuals/setup.md) for requirements. 

## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh kraken2 > SHORT/short_pipeline.out 2>&1 &
```

**Input:**
 - `evaluation/<SampleID>_final_assembly.fasta`
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh kraken2 > LONG/long_pipeline.out 2>&1 &
```

**Input:**
 - `evaluation/<SampleID>_final_assembly.fasta`
 - `clean_reads/<SampleID>_ont.fastq`


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh kraken2 > HYBRID/hybrid_pipeline.out 2>&1 &
```

**Input:**
 - `evaluation/<SampleID>_final_assembly.fasta`
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`
 - `clean_reads/<SampleID>_ont.fastq`

## Output

- `kraken2/<SampleID>/*.krona` - files summarize taxonomy statistics to be fed into KronaTools
- `kraken2/<SampleID>/*.kraken2` - files contain the Kraken-estimated taxonomy of each read or contig
- `kraken2/<SampleID>/kronagram.html` - file contains all the taxonomy information from all the sample and assembly



