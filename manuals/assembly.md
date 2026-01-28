# Metagenomic assembly

**Description:** Assemble pre-processed shotgun metagenomic data, using metaSPAdes (short-reads), metaFlye (long-reads), and OPERA-MS (hybrid)! 


## Module overview

For short-read sequences, this module performs the following tasks:
 1. Assembly using metaSPades
 2. Sorts out unassembled reads
 3. Assembles unassembled reads using Megahit
 4. Combines and formats the two assemblies
 5. Performs assembly QC with QUAST

For long-read sequences, this module performs the following tasks:
 1. Assembly using metaFlye
 2. Performs assembly QC with QUAST

For hybrid assemblies, this module performs the following tasks:
 1. Assembly using OPERA-MS
 2. Performs assembly QC with QUAST

## Setup

See [setup](manuals/setup.md) for requirements. 


## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh assembly > SHORT/short_pipeline.out 2>&1 &
```

**Input:**
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh assembly > LONG/long_pipeline.out 2>&1 &
```

**Input:**
 - `clean_reads/<SampleID>_ont.fastq`


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh assembly > HYBRID/hybrid_pipeline.out 2>&1 &
```

**Input:**
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`
 - `clean_reads/<SampleID>_ont.fastq`


## Output Overview

- `assembly/<SampleID>/assembly.fasta` - assembled metagenome
- `assembly/<SampleID>/assembly_report.html` - the QUAST assembly report
- `assembly/<SampleID>/` - various folders containing intermediate files from this module

