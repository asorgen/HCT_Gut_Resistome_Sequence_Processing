# Sequence trimming and quality filtering

**Description:** Adapter removal and quality trimming of deduplicated Illumina reads using TrimGalore (a wrapper around Cutadapt). Low-quality bases and residual adapter sequences are removed to improve the accuracy of downstream assembly, alignment, and taxonomic classification. A post-trimming FastQC report is generated to verify quality improvement.

## Module overview

### Short-read pipeline

1. Trim adapters and low-quality bases using TrimGalore
2. Rename output files to standardized naming convention
3. Run FastQC quality assessment on trimmed reads

### TrimGalore parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--paired` | | Paired-end mode |
| `--gzip` | | Compress output files |
| `--quality` | 30 | Minimum Phred quality score |
| `--length` | 50 | Minimum read length after trimming |
| `--stringency` | 3 | Minimum overlap with adapter for trimming |
| `--max_n` | 0 | Maximum number of N bases allowed |

## Setup

See [setup](setup.md) for general requirements.

### Software Dependencies

| Software | Source | Notes |
|----------|--------|-------|
| `metawrap-env` conda environment | `anaconda3/2023.09` module | Provides TrimGalore, Cutadapt, and FastQC |

## Configuration

Key variables set in the `.config` files (`pipelineScripts/configs/`):

| Variable | Description | Example |
|----------|-------------|---------|
| `run_trim` | Toggle trimming module on/off | `true` |
| `dedup_Dir` | Input directory (deduplicated reads) | `${datasetDir}/0.2_deduplication` |
| `trimmed_Dir` | Output directory | `${datasetDir}/0.3_sequence_trim` |
| `trim_opts` | SLURM resource options | See table below |

### SLURM resources (`trim_opts`)

| Pipeline | CPUs | Memory | Time |
|----------|------|--------|------|
| Short-read (Duke) | 8 (`--ntasks-per-node`) | 32 GB | 02:00:00 |
| Short-read (UNC) | 8 (`--ntasks-per-node`) | 32 GB | 05:00:00 |

All jobs run on the `Orion` partition with 1 node.

## Run module

The module is launched automatically by the pipeline wrapper when `run_trim=true` in the config. To run the pipeline up through this module, pass the `trim` stop-point argument:

### Short-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh trim > Duke_short_pipeline.out 2>&1 &
```

**Input:**
- `0.2_deduplication/<SampleID>_deduped_R1.fastq.gz`
- `0.2_deduplication/<SampleID>_deduped_R2.fastq.gz`

## Output

### Trimmed reads (`0.3_sequence_trim/`)

| File | Description |
|------|-------------|
| `<SampleID>_trimmed_1.fastq.gz` | Trimmed forward reads |
| `<SampleID>_trimmed_2.fastq.gz` | Trimmed reverse reads |

### FastQC reports (`0.3_sequence_trim/QC_report/`)

| File | Description |
|------|-------------|
| `<SampleID>_trimmed_1_fastqc.html` | FastQC visual report (forward reads) |
| `<SampleID>_trimmed_1_fastqc.zip` | FastQC data archive (forward reads) |
| `<SampleID>_trimmed_2_fastqc.html` | FastQC visual report (reverse reads) |
| `<SampleID>_trimmed_2_fastqc.zip` | FastQC data archive (reverse reads) |

### Completion marker

`0.3_sequence_trim/COMPLETE/<SampleID>` -- created on successful completion.

## Dependencies

**Upstream modules:**
- **0.2 Deduplication** -- provides deduplicated reads from `0.2_deduplication/`

**Downstream modules:**
- **0.4 Host Decontamination** -- uses trimmed reads from `0.3_sequence_trim/`
