# Metagenomic binning

**Description:** This module bins assemblies into metagenome-assembled genomes (MAGs) using three binning algorithms:
 1. MetaBat2
 2. MaxBin2
 3. CONCOCT


## Module overview
 1. Align reads to generate coverage files
  - Index assembly
 2. Run MetaBat2
  - Create contig depth file
  - Perfom binning
 3. Run MaxBin2
  - Create contig depth file
  - Split main contig depth file into inidividual files for input
  - Perform binning
4. Run CONCOCT
 - Index .bam alignment files
 - Cut up contigs
 - Estimate fragment coverage
 - Perform binning
 - Merge 10 kb fragments back into contigs
 - Split contigs into bins

## Setup

See [setup](manuals/setup.md) for requirements. 
 

## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/scripts/short_pipeline.sh binning > SHORT/short_pipeline.out 2>&1 &
```

**Input:**
 - `EVALUATION/<SampleID>_final_assembly.fasta`
 - `CLEAN_READS/<SampleID>_1.fastq`
 - `CLEAN_READS/<SampleID>_2.fastq`


### Long-read analysis

Run command
```shell
nohup sh ./LONG/scripts/long_pipeline.sh binning > LONG/long_pipeline.out 2>&1 &
```

**Input:**
 - `EVALUATION/<SampleID>_final_assembly.fasta`
 - `CLEAN_READS/<SampleID>_ont.fastq`


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/scripts/hybrid_pipeline.sh binning > HYBRID/hybrid_pipeline.out 2>&1 &
```

**Input:**
 - `EVALUATION/<SampleID>_final_assembly.fasta`
 - `CLEAN_READS/<SampleID>_1.fastq`
 - `CLEAN_READS/<SampleID>_2.fastq`
 - `CLEAN_READS/<SampleID>_ont.fastq`


## Output Overview

- `INITIAL_BINNING/concoct_bins` - bins generated using the CONCOCT algorithm
- `INITIAL_BINNING/maxbin2_bins` - bins generated using the MaxBin2 algorithm
- `INITIAL_BINNING/metabat2_bins` - bins generated using the MetaBat2 algorithm
