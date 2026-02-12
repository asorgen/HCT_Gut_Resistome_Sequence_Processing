# Metagenomic Sequence Processing Pipeline

A modular metagenomic analysis pipeline supporting short-read, long-read, and hybrid assembly workflows. The pipeline is orchestrated through SLURM job scheduling with automatic dependency management and uses Conda/metaWRAP environments for reproducibility.

---

## Pipeline Workflow Overview

The pipeline is organized into numbered modules that execute sequentially through SLURM job dependencies. Each module can be toggled on or off via configuration files.

### Preprocessing (Module 0.x)

| Module | Description | Tools | Input | Output |
|--------|-------------|-------|-------|--------|
| 0.1 Pre-QC | Quality assessment and decompression of raw reads | FastQC | Raw gzipped FASTQ files | FastQC reports, unzipped FASTQ files |
| 0.2 Deduplication | Removal of duplicate reads | BBMap (`clumpify.sh`) | Unzipped FASTQ files | De-duplicated FASTQ files |
| 0.3 Trimming | Adapter removal and quality trimming | TrimGalore, Cutadapt | De-duplicated FASTQ files | Quality-trimmed FASTQ files |
| 0.4 Host Decontamination | Removal of human-derived reads | BMTagger (hg38 index) | Trimmed FASTQ files | Clean, host-free FASTQ files |

### Assembly and Evaluation (Module 1.x)

| Module | Description | Tools | Input | Output |
|--------|-------------|-------|-------|--------|
| 1.1 Assembly | *De novo* metagenomic assembly | metaSPAdes (short), metaFlye (long), OPERA-MS (hybrid); Megahit for unassembled reads | Clean FASTQ files | `assembly.fasta`, QUAST report |
| 1.2 Evaluation | Assembly quality filtering | QUAST, custom filtering | `assembly.fasta` | Filtered assembly (contigs >= 1.5 kb), summary statistics |

### Taxonomic Classification (Module 2.x)

| Module | Description | Tools | Input | Output |
|--------|-------------|-------|-------|--------|
| 2.1 Kraken2 | Taxonomic classification of reads and assembly | Kraken2, Bracken, KronaTools | Clean reads + assembly | Kraken2 reports, Bracken abundance estimates, Kronagram HTML visualizations |
| 2.3 MetaPhlAn4 | Marker-based taxonomic profiling (optional) | MetaPhlAn4 | Clean reads | Taxonomic abundance profiles |

### Binning and Refinement (Module 3.x)

| Module | Description | Tools | Input | Output |
|--------|-------------|-------|-------|--------|
| 3.1 Binning | Metagenomic binning using three algorithms | MetaBat2, MaxBin2, CONCOCT | Filtered assembly + clean reads | Three sets of bins (`metabat2_bins/`, `maxbin2_bins/`, `concoct_bins/`) |
| 3.2 Refinement | Bin consolidation and quality control | metaWRAP `refine_bins`, CheckM | All binning outputs | Refined bins passing quality thresholds (default: >= 70% completion, <= 10% contamination) |
| 3.3 Reassembly | Re-assembly of bin-associated reads | SPAdes | Refined bins + clean reads | Reassembled MAGs (metagenome-assembled genomes) |

### Functional Annotation (Module 4.x)

| Module | Description | Tools | Input | Output |
|--------|-------------|-------|-------|--------|
| 4.1 Classify Bins | Taxonomic assignment of MAGs | BLASTN, NCBI databases | Reassembled bins | Bin-level taxonomic classifications |
| 4.2 Annotate Bins | Gene prediction and functional annotation | Prokka, Bakta | Reassembled bins | GFF annotation files, translated/untranslated gene sequences |

### AMR Detection (Module 5.x)

Multiple complementary approaches for antimicrobial resistance gene identification:

| Module | Description | Tools | Input | Output |
|--------|-------------|-------|-------|--------|
| 5.1 NT AMR (Assembly) | Nucleotide-level AMR detection from assembly | AMRFinder+, RGI | Filtered assembly | AMR gene predictions |
| 5.2 NT AMR (MAGs) | Nucleotide-level AMR detection from MAGs | AMRFinder+, RGI | Reassembled bins | Bin-level AMR predictions |
| 5.3 ShortBRED | Marker-based AMR quantification from reads | ShortBRED (CARD markers) | Clean reads | AMR marker abundances |
| 5.4 RGI BWT | Read-mapping AMR profiling against CARD | RGI (KMA aligner) | Clean reads | CARD-mapped AMR profiles |
| 5.5 AA AMR (Assembly) | Protein-level AMR detection from assembly | Prodigal, AMRFinder+, RGI | Filtered assembly | Protein-level AMR predictions |
| 5.6 AA AMR (MAGs) | Protein-level AMR detection from MAGs | Prodigal, AMRFinder+, RGI | Reassembled bins | Protein-level bin AMR predictions |

### Cleanup (Module 6)

Removes intermediate files after all modules complete to reclaim storage.

---

## Directory Structure

```
sequence_processing/
├── README.md
├── manuals/                            # Detailed documentation for each module
├── pipelineScripts/
│   ├── configs/                        # Pipeline configuration files (.config)
│   │   └── functions.sh                # Shared helper functions (job submission, logging)
│   ├── pipeline_wrappers/              # Entry-point scripts that orchestrate the full pipeline
│   ├── sample_lists/                   # Sample ID manifests (one sample per line)
│   └── scripts/
│       ├── short_scripts/              # Short-read pipeline modules (0.x - 6)
│       ├── long_scripts/               # Long-read pipeline modules
│       ├── hybrid_scripts/             # Hybrid assembly modules
│       ├── helper_scripts/             # Utility scripts (parsing, counting, database download)
│       └── post_scripts/               # Post-processing summarization and aggregation
└── databases/                          # Reference databases (Kraken2, CARD, CheckM, etc.)
```

---

## Configuration

Each pipeline variant is controlled by a `.config` file in `pipelineScripts/configs/`. Configuration files define:

- **Data paths**: Locations for raw sequences, processed output, and databases
- **Module toggles**: Enable or disable individual pipeline modules (YES/NO flags)
- **SLURM resources**: Cores, memory, and walltime for each module
- **Quality thresholds**: Bin completion/contamination cutoffs, minimum contig length
- **Database paths**: Paths to Kraken2, CheckM, Bakta, CARD, and other reference databases

Available configuration files:

| Config File | Description |
|-------------|-------------|
| `Duke_short-read.config` | Short-read (Illumina) pipeline |
| `Duke_long-read.config` | Long-read (Oxford Nanopore) pipeline |
| `Duke_hybrid-read.config` | Hybrid (short + long read) assembly pipeline |
| `UNC_short-read.config` | Short-read pipeline (alternate cohort) |

---

## Pipeline Modes

The pipeline supports two execution modes, selectable in the configuration file:

### Full Mode
Runs the complete pipeline from preprocessing through protein-level AMR detection, including assembly, binning, bin refinement, reassembly, functional annotation, and all AMR detection approaches.

### Heston Mode
A streamlined mode that runs only preprocessing and read-based AMR profiling. This skips assembly, binning, and annotation modules, enabling faster turnaround for resistance gene profiling:
- **Enabled**: Pre-QC, Deduplication, Trimming, Host Decontamination, Kraken2, RGI BWT
- **Disabled**: Assembly, Evaluation, Binning, Refinement, Reassembly, Classification, Annotation, ShortBRED, assembly/MAG-based AMR

---

## Running the Pipeline

### Prerequisites
1. Access to a SLURM workload manager
2. Conda environments with required tools installed (see [Setup](manuals/setup.md))
3. Required reference databases downloaded to the `databases/` directory
4. A sample list file with one sample ID per line in `pipelineScripts/sample_lists/`
5. A completed `.config` file with correct paths and module settings

### Quickstart

Pipeline wrapper scripts are located in `pipelineScripts/pipeline_wrappers/`. Each wrapper reads its corresponding config file, iterates through samples, and submits SLURM jobs with dependency chains.

```bash
# Short-read pipeline
nohup bash pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh > short_pipeline.out 2>&1 &

# Long-read pipeline
nohup bash pipelineScripts/pipeline_wrappers/Duke_long_pipeline.sh > long_pipeline.out 2>&1 &

# Hybrid assembly pipeline
nohup bash pipelineScripts/pipeline_wrappers/Duke_hybrid_pipeline.sh > hybrid_pipeline.out 2>&1 &
```

The wrapper scripts handle:
- Reading the configuration file
- Iterating through all samples in the sample list
- Submitting each module as a SLURM job with appropriate dependencies
- Checking for already-completed modules (via flag files) and skipping them
- Logging job IDs and submission status

---

## Module Manuals

Detailed documentation for each pipeline module:

1. [Setup](manuals/setup.md)
2. [**Pre-processing** metagenomic data](manuals/preprocessing.md)
3. [Metagenomic **assembly**](manuals/assembly.md)
4. [Assembly refinement and **evaluation**](manuals/evaluation.md)
5. [Taxonomic classification with **Kraken2**](manuals/kraken2.md)
6. [Metagenomic **binning**](manuals/binning.md)
7. [Metagenomic bin **refinement**](manuals/refine.md)
8. [Metagenomic bin **reassembly**](manuals/reassembly.md)
9. [Taxonomic bin **classification**](manuals/classify_bins.md)
10. [**Functional annotation** of bins](manuals/annotation.md)
11. [**AMR detection**](manuals/amr_search.md)

---

## Key Tools and Dependencies

| Category | Tools |
|----------|-------|
| Quality Control | FastQC, QUAST |
| Read Processing | BBMap, TrimGalore, Cutadapt, BMTagger |
| Assembly | metaSPAdes, Megahit, metaFlye, OPERA-MS, SPAdes |
| Taxonomic Classification | Kraken2, Bracken, KronaTools, MetaPhlAn4, BLASTN |
| Binning | MetaBat2, MaxBin2, CONCOCT |
| Bin Refinement | metaWRAP, CheckM |
| Annotation | Prokka, Bakta, Prodigal |
| AMR Detection | AMRFinder+, RGI (CARD), ShortBRED |

---

## Required Databases

| Database | Purpose |
|----------|---------|
| hg38 (BMTagger index) | Human host read removal |
| Kraken2 DB | Taxonomic classification |
| CARD | Comprehensive Antibiotic Resistance Database (used by RGI, ShortBRED) |
| CheckM | Bin quality assessment |
| GTDBtk | Genome taxonomy classification |
| Bakta | Functional annotation |
| ShortBRED markers | CARD-derived AMR markers for read-based quantification |
| NCBI BLAST DB | Nucleotide similarity searches |
| MetaPhlAn4 markers | Marker-based taxonomic profiling |
