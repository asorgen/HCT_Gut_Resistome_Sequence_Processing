# Host decontamination

**Description:** Removal of human-derived reads from metagenomic samples by aligning against the GRCh38 human reference genome with Bowtie2. Reads that map to the human genome are discarded, retaining only microbial reads for downstream analysis. A post-decontamination FastQC report is generated to assess the final cleaned read set. Intermediate alignment files are automatically cleaned up after completion.

## Module overview

### Short-read pipeline

1. Index the GRCh38 human reference genome with Bowtie2 (first run only)
2. Decompress trimmed reads
3. Align reads to the human genome (Bowtie2)
4. Convert SAM to BAM format (samtools view)
5. Filter for unmapped reads only (samtools view `-f 12 -F 256`)
6. Sort BAM by read name for proper pairing (samtools sort)
7. Extract paired reads to FASTQ files (samtools fastq)
8. Run FastQC on final cleaned reads
9. Remove intermediate files (unzipped reads, SAM, BAM)

## Setup

See [setup](setup.md) for general requirements.

### Software Dependencies

| Software | Source | Notes |
|----------|--------|-------|
| `metawrap-env` conda environment | `anaconda3/2023.09` module | Provides FastQC |
| Bowtie2 | `module load bowtie2` | Read alignment |
| Samtools | `module load samtools` | SAM/BAM processing |

### Database Requirements

| Database | Path | Description |
|----------|------|-------------|
| GRCh38 reference genome | `${GRCh38_DB}/GRCh38.p14.fa` | Human genome reference (patch 14) |
| Bowtie2 index | `${GRCh38_DB}/bowtie2_index/GRCh38.p14` | Built automatically on first run |

## Configuration

Key variables set in the `.config` files (`pipelineScripts/configs/`):

| Variable | Description | Example |
|----------|-------------|---------|
| `run_decontam` | Toggle decontamination module on/off | `true` |
| `trimmed_Dir` | Input directory (trimmed reads) | `${datasetDir}/0.3_sequence_trim` |
| `clean_readDir` | Output directory | `${datasetDir}/0.4_host_decontamination` |
| `GRCh38_DB` | Path to GRCh38 reference genome directory | `${ROOT}/databases/GRCh38` |
| `decontam_opts` | SLURM resource options | See table below |

### SLURM resources (`decontam_opts`)

| Pipeline | CPUs | Memory | Time |
|----------|------|--------|------|
| Short-read (Duke) | 32 (`--ntasks-per-node`) | 136 GB | 02:00:00 |
| Short-read (UNC) | 32 (`--ntasks-per-node`) | 136 GB | 02:00:00 |

All jobs run on the `Orion` partition with 1 node.

## Run module

The module is launched automatically by the pipeline wrapper when `run_decontam=true` in the config. To run the pipeline up through this module, pass the `decontam` stop-point argument:

### Short-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh decontam > Duke_short_pipeline.out 2>&1 &
```

**Input:**
- `0.3_sequence_trim/<SampleID>_trimmed_1.fastq.gz`
- `0.3_sequence_trim/<SampleID>_trimmed_2.fastq.gz`
- `${GRCh38_DB}/GRCh38.p14.fa` (human reference genome)

## Output

### Cleaned reads (`0.4_host_decontamination/`)

| File | Description |
|------|-------------|
| `<SampleID>_1.fastq.gz` | Final cleaned forward reads (host-free) |
| `<SampleID>_2.fastq.gz` | Final cleaned reverse reads (host-free) |

### FastQC reports (`0.4_host_decontamination/QC_report/`)

| File | Description |
|------|-------------|
| `<SampleID>_1_fastqc.html` | FastQC visual report (forward reads) |
| `<SampleID>_1_fastqc.zip` | FastQC data archive (forward reads) |
| `<SampleID>_2_fastqc.html` | FastQC visual report (reverse reads) |
| `<SampleID>_2_fastqc.zip` | FastQC data archive (reverse reads) |

### Completion marker

`0.4_host_decontamination/COMPLETE/<SampleID>` -- created on successful completion. Intermediate files (unzipped reads, SAM, BAM, sorted BAM) are removed at this point.

## Dependencies

**Upstream modules:**
- **0.3 Sequence Trimming** -- provides trimmed reads from `0.3_sequence_trim/`

**Downstream modules:**
- **1.1 Assembly** -- uses cleaned reads from `0.4_host_decontamination/`
- **2.1 Kraken2** -- uses cleaned reads for taxonomic classification
- **5.3 ShortBRED** -- uses cleaned reads for AMR marker quantification
- **5.4 RGI BWT** -- uses cleaned reads for read-mapping AMR profiling
