# Taxonomic classification with Kraken2

**Description**: Running Kraken2 on the reads will give us an idea of the taxonomic composition of the community. Running Kraken2 on the assembly will give us an idea what taxonomic groups were assembled better than others (the assembly process is heavily biased and should not be used to infer overall community composition). Bracken then re-estimates species-level abundances from the Kraken2 reports, correcting for unclassified reads at higher taxonomic levels.

## Module overview

### Short-read pipeline

1. Unzip cleaned sequence files (`.fastq.gz` to `.fastq`)
2. Run Kraken2 on reads and assembly via the modified metaWRAP `kraken2_bracken` module
3. Translate Kraken2 output to taxonomy format
4. Generate kronagrams with KronaTools
5. Estimate species-level abundances with Bracken

### Long-read / Hybrid pipeline

1. Run Kraken2 directly (no metaWRAP wrapper) on ONT reads, assembly, and short reads (hybrid only)
2. Translate Kraken2 output to taxonomy format
3. Generate kronagrams with KronaTools
4. Estimate species-level abundances with Bracken

## Setup

See [setup](setup.md) for general requirements.

### Modified metaWRAP Module Installation

This module uses a modified version of metaWRAP's `kraken2` module called `kraken2_bracken`, which integrates Kraken2 classification, KronaTools visualization, and Bracken abundance estimation into a single workflow. The modified module is included in this repository under `pipelineScripts/modified-metaWRAP-modules/`.

To install it into your metaWRAP installation:

1. Locate your metaWRAP installation's `bin/` directory (where the `metawrap` script lives):
   ```bash
   which metawrap
   # Example output: /path/to/metaWRAP/bin/metawrap
   ```

2. Copy the modified module into metaWRAP's `bin/metawrap-scripts/` directory:
   ```bash
   METAWRAP_DIR="$(dirname "$(which metawrap)")/.."
   cp pipelineScripts/modified-metaWRAP-modules/kraken2_bracken.sh "$METAWRAP_DIR/bin/metawrap-scripts/"
   ```

3. Register the module by adding the following block to the `metawrap` master script (located at `$METAWRAP_DIR/bin/metawrap`). Add it after the existing `kraken2` entry:
   ```bash
   elif [ "$1" = kraken2_bracken ]; then
   	echo metawrap kraken2_bracken ${@:2}
   	time ${PIPES}/kraken2_bracken.sh ${@:2}
   ```

4. Verify the installation:
   ```bash
   metawrap kraken2_bracken --help
   ```

### Database Requirements

| Database | Path | Description |
|----------|------|-------------|
| Kraken2 DB | `${ROOT}/databases/KRAKEN2-TESSA-DB` | TESSA Kraken2 reference database (~125 GB) |
| Bracken k-mer distribution | `${KRAKEN2_DB}/database150mers.kmer_distrib` | Pre-built 150-mer distribution file inside the Kraken2 DB |

### Software Dependencies

| Software | Source | Used by |
|----------|--------|---------|
| `metawrap-env` conda environment | `anaconda3/2023.09` module | Short-read pipeline |
| Bracken 2.7 (`est_abundance.py`) | `${bracken}` config variable | All pipelines |
| KronaTools (`ktImportText`) | Included in metaWRAP env | All pipelines |

## Configuration

Key variables set in the `.config` files (`pipelineScripts/configs/`):

| Variable | Description | Example |
|----------|-------------|---------|
| `run_k2` | Toggle Kraken2 module on/off | `true` |
| `KRAKEN2_DB` | Path to Kraken2 database | `${ROOT}/databases/KRAKEN2-TESSA-DB` |
| `bracken` | Path to Bracken `src/` directory | `/users/asorgen/PROGRAMS/Bracken-2.7/src` |
| `krakenDir` | Kraken2 output directory | `${datasetDir}/2.1_kraken2` |
| `brackenDir` | Bracken output directory | `${datasetDir}/2.2_bracken` |
| `k2_opts` | SLURM resource options for this module | See table below |

### SLURM resources (`k2_opts`)

| Pipeline | CPUs | Memory | Time | Notes |
|----------|------|--------|------|-------|
| Short-read | 16 (`--cpus-per-task`) | 256 GB | 00:45:00 | |
| Long-read | 16 (`--ntasks-per-node`) | 256 GB | 01:00:00 | |
| Hybrid | 16 (`--ntasks-per-node`) | 256 GB | 24:00:00 | Processes short reads, ONT reads, and assembly |

All jobs run on the `Orion` partition with 1 node.

## Run module

The module is launched automatically by the pipeline wrapper when `run_k2=true` in the config. To run the pipeline up through this module, pass the `kraken2` stop-point argument:

### Short-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh kraken2 > Duke_short_pipeline.out 2>&1 &
```

**Input:**
 - `clean_reads/<SampleID>_1.fastq.gz`
 - `clean_reads/<SampleID>_2.fastq.gz`
 - `evaluation/<SampleID>_final_assembly.fasta` (optional)

### Long-read analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_long_pipeline.sh kraken2 > Duke_long_pipeline.out 2>&1 &
```

**Input:**
 - `clean_reads/<SampleID>_ont.fastq`
 - `evaluation/<SampleID>_final_assembly.fasta`

### Hybrid assembly analysis

```shell
nohup sh ./pipelineScripts/pipeline_wrappers/Duke_hybrid_pipeline.sh kraken2 > Duke_hybrid_pipeline.out 2>&1 &
```

**Input:**
 - `clean_reads/<SampleID>_1.fastq.gz`
 - `clean_reads/<SampleID>_2.fastq.gz`
 - `clean_reads/<SampleID>_ont.fastq`
 - `evaluation/<SampleID>_final_assembly.fasta`

## Output

### Kraken2 output (`2.1_kraken2/<SampleID>/`)

| File | Description |
|------|-------------|
| `<SampleID>.krak2` | Raw Kraken2 classification output (reads) |
| `<SampleID>.kreport` | Kraken2 report format (reads) |
| `<SampleID>.kraken2` | Translated taxonomy |
| `<SampleID>.krona` | KronaTools input format |
| `<SampleID>_assembly.krak2` | Raw Kraken2 classification output (assembly) |
| `<SampleID>_assembly.kreport` | Kraken2 report format (assembly) |
| `<SampleID>_assembly.kraken2` | Translated taxonomy (assembly) |
| `<SampleID>_assembly.krona` | KronaTools input format (assembly) |
| `kronagram.html` | Combined interactive Krona visualization |

Long-read and hybrid pipelines also produce `<SampleID>_ONT.*` files for ONT reads.

### Bracken output (`2.2_bracken/`)

| File | Description |
|------|-------------|
| `sr/<SampleID>.bracken.out` | Species-level abundance from short reads |
| `assembly/<SampleID>_assembly.bracken.out` | Species-level abundance from assembly |
| `ont/<SampleID>_ONT.bracken.out` | Species-level abundance from ONT reads (long-read/hybrid only) |

### Completion marker

`2.1_kraken2/COMPLETE/<SampleID>` — created on successful completion. Intermediate files (unzipped fastqs) are cleaned up at this point.

## Dependencies

**Upstream modules:**
- **0.4 Host Decontamination** — provides cleaned reads (`clean_reads/`)
- **1.2 Evaluation** — provides the final assembly (`evaluation/<SampleID>_final_assembly.fasta`); optional for short-read Kraken2 on reads, required for assembly classification
