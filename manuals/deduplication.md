# Sequence deduplication

**Description:** Deduplication removes PCR and optical duplicate reads from Illumina sequencing data using HTStream's SuperDeduper. Duplicates inflate apparent coverage and can bias downstream analyses including assembly, binning, and abundance estimation. A post-deduplication FastQC report is generated to assess the effect of duplicate removal.

## Module overview

### Short-read pipeline

1. Remove duplicate reads using HTStream SuperDeduper
2. Run FastQC quality assessment on deduplicated reads

## Setup

See [setup](setup.md) for general requirements.

### Software Dependencies

| Software | Source | Notes |
|----------|--------|-------|
| `htstream-env` conda environment | `anaconda3/2023.09` module | Provides `hts_SuperDeduper` |
| `metawrap-env` conda environment | `anaconda3/2023.09` module | Provides FastQC |

## Configuration

Key variables set in the `.config` files (`pipelineScripts/configs/`):

| Variable | Description | Example |
|----------|-------------|---------|
| `run_dedup` | Toggle deduplication module on/off | `true` |
| `raw_readDir` | Input directory (raw reads) | `${datasetDir}/0.0_raw_reads` |
| `dedup_Dir` | Output directory | `${datasetDir}/0.2_deduplication` |
| `dedup_opts` | SLURM resource options | See table below |

### SLURM resources (`dedup_opts`)

| Pipeline | CPUs | Memory | Time |
|----------|------|--------|------|
| Short-read (Duke) | 16 (`--ntasks-per-node`) | 64 GB | 01:30:00 |
| Short-read (UNC) | 16 (`--ntasks-per-node`) | 64 GB | 01:30:00 |

All jobs run on the `Orion` partition with 1 node.

## Run module

The module is launched automatically by the pipeline wrapper when `run_dedup=true` in the config. To run the pipeline up through this module, pass the `dedup` stop-point argument:

### Short-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh dedup > Duke_short_pipeline.out 2>&1 &
```

**Input:**
- `0.0_raw_reads/<SampleID>_1.fastq`
- `0.0_raw_reads/<SampleID>_2.fastq`

## Output

### Deduplication output (`0.2_deduplication/`)

| File | Description |
|------|-------------|
| `<SampleID>_deduped_R1.fastq.gz` | Deduplicated forward reads |
| `<SampleID>_deduped_R2.fastq.gz` | Deduplicated reverse reads |
| `<SampleID>_stats.log` | SuperDeduper statistics (duplicate counts, rates) |

### FastQC reports (`0.2_deduplication/QC_report/`)

| File | Description |
|------|-------------|
| `<SampleID>_deduped_R1_fastqc.html` | FastQC visual report (forward reads) |
| `<SampleID>_deduped_R1_fastqc.zip` | FastQC data archive (forward reads) |
| `<SampleID>_deduped_R2_fastqc.html` | FastQC visual report (reverse reads) |
| `<SampleID>_deduped_R2_fastqc.zip` | FastQC data archive (reverse reads) |

### Completion marker

`0.2_deduplication/COMPLETE/<SampleID>` -- created on successful completion.

## Dependencies

**Upstream modules:**
- **0.1 Pre-QC** -- provides raw reads in `0.0_raw_reads/`

**Downstream modules:**
- **0.3 Sequence Trimming** -- uses deduplicated reads from `0.2_deduplication/`
