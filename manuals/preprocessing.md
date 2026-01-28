# Pre-processing metagenomic data

**Description:** Pre-processing raw metagenomic data is necessary before any other application. 


## Module overview

For short-read sequences, this module performs the following tasks:
 1. Run pre-QC report (FastQC)
 2. Adapter trimming (cutadapt) + quality trimming with TrimGalore
 3. Remove host sequence contamination (bmtagger)
 4. Run post-QC report (FastQC)

For long-read sequences, this module performs the following tasks:
 1. Run pre-QC report (FastQC)
 2. Adapter trimming (cutadapt) + quality trimming with Porechop
 3. Remove host sequence contamination (bmtagger)
 4. Run post-QC report (FastQC)


## Setup

See [setup](manuals/setup.md) for requirements. 


## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh read_qc > SHORT/short_pipeline.out 2>&1 &
```

**Input:**
 - raw forward reads
 - raw reverse reads

**Output:**
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`



### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh read_qc > LONG/long_pipeline.out 2>&1 &
```

**Input:**
 - raw Nanopore reads

**Output:**
 - `clean_reads/<SampleID>_ont.fastq`


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh read_qc > HYBRID/hybrid_pipeline.out 2>&1 &
```

**Input:**
 - raw forward reads
 - raw reverse reads
 - raw Nanopore reads

**Output:**
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`
 - `clean_reads/<SampleID>_ont.fastq`

## Output Overview

**This module creates the following output in your SHORT, LONG, OR HYBRID folders:**
- `raw_reads/` - a temporary storage place for your raw sequence files
- `read_qc/<SampleID>/pre-QC_report`- quality report of raw sequences
- `read_qc/<SampleID>/post-QC_report` - quality report of clean, decontaminated sequences
- `clean_reads/` location of quality-filtered, decontaminated sequence files



