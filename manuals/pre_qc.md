# Pre-QC quality assessment

**Description:** Pre-QC performs quality assessment and decompression of raw metagenomic sequence reads. Raw gzipped FASTQ files are decompressed, renamed with standardized sample IDs, and assessed with FastQC to evaluate base quality, adapter content, duplication rates, and other quality metrics before downstream processing.

## Module overview

### Short-read pipeline

1. Decompress and rename raw sequence files (`gunzip`)
2. Run FastQC quality assessment on raw reads

## Setup

See [setup](setup.md) for general requirements.

### Software Dependencies

| Software | Source | Notes |
|----------|--------|-------|
| `metawrap-env` conda environment | `anaconda3/2023.09` module | Provides FastQC |

## Configuration

Key variables set in the `.config` files (`pipelineScripts/configs/`):

| Variable | Description | Example |
|----------|-------------|---------|
| `run_pre_qc` | Toggle pre-QC module on/off | `true` |
| `seqPath` | Path to raw gzipped sequence files | See config |
| `R1_ext` | Forward read file extension pattern | `R1_001.fastq.gz` |
| `R2_ext` | Reverse read file extension pattern | `R2_001.fastq.gz` |
| `raw_readDir` | Directory for decompressed raw reads | `${datasetDir}/0.0_raw_reads` |
| `pre_qcDir` | Pre-QC output directory | `${datasetDir}/0.1_pre_qc` |
| `pre_qc_opts` | SLURM resource options | See table below |

### SLURM resources (`pre_qc_opts`)

| Pipeline | CPUs | Memory | Time |
|----------|------|--------|------|
| Short-read (Duke) | 8 (`--ntasks-per-node`) | 16 GB | 00:30:00 |
| Short-read (UNC) | 8 (`--ntasks-per-node`) | 16 GB | 00:30:00 |

All jobs run on the `Orion` partition with 1 node.

## Run module

The module is launched automatically by the pipeline wrapper when `run_pre_qc=true` in the config. To run the pipeline up through this module, pass the `pre_qc` stop-point argument:

### Short-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh pre_qc > Duke_short_pipeline.out 2>&1 &
```

**Input:**
- `${seqPath}/<IlluminaID>_${R1_ext}` (raw gzipped forward reads)
- `${seqPath}/<IlluminaID>_${R2_ext}` (raw gzipped reverse reads)

## Output

### Pre-QC output (`0.1_pre_qc/`)

| File | Description |
|------|-------------|
| `<SampleID>_1_fastqc.html` | FastQC visual report (forward reads) |
| `<SampleID>_1_fastqc.zip` | FastQC data archive (forward reads) |
| `<SampleID>_2_fastqc.html` | FastQC visual report (reverse reads) |
| `<SampleID>_2_fastqc.zip` | FastQC data archive (reverse reads) |

### Raw reads (`0.0_raw_reads/`)

| File | Description |
|------|-------------|
| `<SampleID>_1.fastq` | Decompressed, renamed forward reads |
| `<SampleID>_2.fastq` | Decompressed, renamed reverse reads |

### Completion marker

`0.1_pre_qc/COMPLETE/<SampleID>` -- created on successful completion.

## Dependencies

**Upstream modules:**
- None -- this is the first module in the pipeline

**Downstream modules:**
- **0.2 Deduplication** -- uses raw reads from `0.0_raw_reads/`
