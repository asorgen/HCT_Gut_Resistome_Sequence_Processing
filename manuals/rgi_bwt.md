# RGI BWT read-mapping AMR profiling

**Description:** RGI BWT (Resistance Gene Identifier - Burrows-Wheeler Transform) performs read-level antimicrobial resistance (AMR) gene detection by aligning quality-filtered metagenomic reads directly against the CARD (Comprehensive Antibiotic Resistance Database) reference sequences using the KMA aligner. Unlike assembly-based approaches, RGI BWT captures AMR genes from reads that may not assemble well, providing FPKM-normalized abundance estimates for resistance genes across samples.

## Module overview

### Short-read pipeline

1. Load the local CARD database (first run per sample only)
   - Import `card.json`
   - Generate CARD annotation FASTA
   - Load annotations into local RGI database
2. Run `rgi bwt` with KMA aligner on paired-end cleaned reads
3. Validate output (check for hits in gene mapping data)
4. Copy gene mapping results to project output directory
5. Clean up temporary scratch directory

## Setup

See [setup](setup.md) for general requirements.

### Software Dependencies

| Software | Source | Notes |
|----------|--------|-------|
| RGI | `${RGI_ENV}` (Python virtualenv) | Resistance Gene Identifier |
| KMA | Bundled with RGI | K-mer alignment for read mapping |
| Diamond 2.0.9 | `module load diamond/2.0.9` | Protein sequence aligner |
| Samtools | `module load samtools` | SAM/BAM processing |
| Bamtools | `module load bamtools` | BAM utilities |
| Bedtools2 | `module load bedtools2` | Genomic interval operations |

### Database Requirements

| Database | Path | Description |
|----------|------|-------------|
| CARD database | `${CARD_DB}/card.json` | CARD reference database (JSON format) |
| CARD DB directory | `${CARD_DB}` | Working directory for RGI on scratch storage |

The CARD database version should match between the sequence processing and downstream analysis pipelines. The current pipeline uses CARD v4.0.1.

## Configuration

Key variables set in the `.config` files (`pipelineScripts/configs/`):

| Variable | Description | Example |
|----------|-------------|---------|
| `run_rgi_bwt` | Toggle RGI BWT module on/off | `true` |
| `clean_readDir` | Input directory (cleaned reads) | `${datasetDir}/0.4_host_decontamination` |
| `rgi_bwt_dir` | Output directory | `${datasetDir}/5.4_rgi_bwt` |
| `RGI_ENV` | Path to RGI Python virtualenv | `${HPC_HOME}/PROGRAMS/rgi_env/bin/activate` |
| `CARD_DB` | Path to CARD database on scratch | `${HPC_SCRATCH}/RGI_databases/rgi_4.0.1` |
| `rgi_bwt_opts` | SLURM resource options | See table below |

### SLURM resources (`rgi_bwt_opts`)

| Pipeline | CPUs | Memory | Time |
|----------|------|--------|------|
| Short-read (Duke) | 16 (`--cpus-per-task`) | 64 GB | 3:00:00 |
| Short-read (UNC) | 8 (`--cpus-per-task`) | 128 GB | 3:00:00 |

All jobs run on the `Orion` partition with 1 node.

## Run module

The module is launched automatically by the pipeline wrapper when `run_rgi_bwt=true` in the config. In Heston mode, this module runs after host decontamination and Kraken2.

### Short-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh > Duke_short_pipeline.out 2>&1 &
```

To run only through RGI BWT, ensure `run_rgi_bwt=true` in the config and set downstream modules to `false`.

**Input:**
- `0.4_host_decontamination/<SampleID>_1.fastq.gz` (cleaned forward reads)
- `0.4_host_decontamination/<SampleID>_2.fastq.gz` (cleaned reverse reads)
- `${CARD_DB}/card.json` (CARD reference database)

## Output

### RGI BWT output (`5.4_rgi_bwt/kma_output/`)

| File | Description |
|------|-------------|
| `<SampleID>.rgi_kma.txt` | Gene mapping data with AMR gene hits, alignment statistics, and abundance metrics |

The gene mapping output file contains per-gene information including:
- ARO (Antibiotic Resistance Ontology) accession
- Gene name and AMR gene family
- Drug class and resistance mechanism
- Read counts and coverage metrics (used downstream for FPKM normalization)

### Completion marker

`5.4_rgi_bwt/COMPLETE/<SampleID>` -- created on successful completion. Temporary per-sample scratch directories (`${CARD_DB}/<SampleID>_kma/`) are cleaned up at this point.

### Temporary files

RGI BWT creates a per-sample working directory on scratch storage (`${CARD_DB}/<SampleID>_kma/`) during execution. This directory contains:
- Local CARD database copy (`localDB/`)
- KMA alignment index and intermediate files
- Raw gene mapping output

These files are automatically removed after the gene mapping data is successfully copied to the project output directory.

## Dependencies

**Upstream modules:**
- **0.4 Host Decontamination** -- provides cleaned reads from `0.4_host_decontamination/`

**Downstream modules:**
- **Post-pipeline summary** (`pipeline_summary.sh`) -- aggregates per-sample RGI BWT output into combined FPKM tables used by the Heston reproduction analysis pipeline
