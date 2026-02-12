# Setup

This guide covers installing all tools and databases needed to run the pipeline from scratch.

---

## Prerequisites

- A SLURM-managed HPC cluster with `module load` support
- Conda (via Anaconda or Miniforge)
- Sufficient storage for reference databases (~200+ GB)

---

## 1. Conda Environments

The pipeline relies on several Conda environments. Create them as follows:

### metaWRAP environment (core pipeline)

metaWRAP is used throughout the pipeline for assembly, binning, refinement, taxonomic classification, and more.

```shell
conda create -n metawrap-env -c ursky -c bioconda -c conda-forge -c defaults \
  metawrap-mg=1.3 \
  fastqc trimgalore cutadapt \
  megahit spades quast \
  metabat2 maxbin2 concoct \
  prokka bakta prodigal \
  kraken2 kronatools \
  blast
```

> **Note:** After installation, locate the metaWRAP config file with `which config-metawrap` and update the database paths as described in the [Databases](#2-databases) section below.

### HTStream environment (deduplication)

```shell
conda create -n htstream-env -c bioconda htstream
```

### Bakta environment (functional annotation)

```shell
conda create -n bakta-env -c bioconda bakta
```

### AMRFinder+ environment

```shell
# Using Miniforge/Mamba is recommended for AMRFinder+
mamba create -n amrfinderplus-4.0.3 -c bioconda -c conda-forge ncbi-amrfinderplus=4.0.3
```

### RGI environment (CARD-based AMR detection)

```shell
# RGI is installed in a virtual environment
python -m venv rgi_env
source rgi_env/bin/activate
pip install rgi
deactivate
```

## HPC Modules

The following tools are expected to be available via `module load` on your cluster. Check availability with `module avail <tool>`:

| Module | Used By |
|--------|---------|
| `anaconda3` | All modules (Conda initialization) |
| `bowtie2` | 0.4 Host Decontamination |
| `samtools` | 0.4 Host Decontamination, 5.4 RGI BWT |
| `multiqc` | Pipeline summary (read QC aggregation) |
| `metaphlan/4.2.2` | 2.3 MetaPhlAn4 |
| `gtdbtk/2.4.0` | 4.1 Classify Bins |
| `blast/2.11.0+` | 5.1/5.2/5.5/5.6 AMR detection |
| `hmmer/3.3.2` | 5.1/5.2/5.5/5.6 AMR detection |
| `diamond/2.0.9` | 5.1/5.2/5.5/5.6 AMR detection, 5.4 RGI BWT |
| `bamtools` | 5.4 RGI BWT |
| `bedtools2` | 5.4 RGI BWT |
| `shortbred/0.9.5` | 5.3 ShortBRED |
| `R/4.3.3` | Post-processing (MetaPhlAn4 tables) |
| `flye` | Long-read assembly |

If any module is not available on your cluster, install it via Conda or from source and update the scripts accordingly.

---

## 2. Databases

All databases are stored in the `databases/` directory under the pipeline root. Paths are configured in the `.config` files.

### GRCh38 Human Reference Genome (Bowtie2 index)

Used by module 0.4 for host read removal via Bowtie2 alignment.

```shell
mkdir -p databases/GRCh38
cd databases/GRCh38

# Download the GRCh38 reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz
gunzip GCA_000001405.29_GRCh38.p14_genomic.fna.gz
mv GCA_000001405.29_GRCh38.p14_genomic.fna GRCh38.p14.fa

# Build the Bowtie2 index (this requires significant memory and time)
mkdir bowtie2_index
bowtie2-build GRCh38.p14.fa bowtie2_index/GRCh38.p14
```

> **Note:** The 0.4 script will automatically build the Bowtie2 index if the `bowtie2_index/` directory does not exist, but pre-building is recommended.

### Kraken2 Taxonomic Database

Any standard Kraken2 database will work. To build a custom database:

```shell
mkdir -p databases/KRAKEN2-DB
kraken2-build --standard --db databases/KRAKEN2-DB --threads 16

# Build Bracken k-mer distribution (required for species-level abundance estimation)
bracken-build -d databases/KRAKEN2-DB -t 16 -k 35 -l 150
```

The Bracken build step creates `database150mers.kmer_distrib` inside the Kraken2 database directory, which is required by module 2.1.

**Don't forget to specify the `KRAKEN2_DB` variable in the metaWRAP config file.** Run `which config-metawrap` to find it.

### CheckM Database

Required by metaWRAP's bin refinement module (3.2).

```shell
mkdir -p databases/MY_CHECKM_FOLDER

# Download CheckM data
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xzf checkm_data_2015_01_16.tar.gz -C databases/MY_CHECKM_FOLDER

# Tell CheckM where to find the data
checkm data setRoot databases/MY_CHECKM_FOLDER
```

### GTDBtk Database

Required by module 4.1 for genome taxonomy classification of MAGs.

```shell
mkdir -p databases/GTDBtk
cd databases/GTDBtk

# Download the latest GTDB-Tk reference data (check https://ecogenomics.github.io/GTDBTk/ for current release)
wget https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar -xzf gtdbtk_data.tar.gz
```

Set the `GTDBTK_DATA_PATH` environment variable or update the config to point to the extracted directory.

### Bakta Database

Required by module 4.2 for functional annotation.

```shell
conda activate bakta-env
bakta_db download --output databases/bakta/db
conda deactivate
```

### CARD Database (RGI)

Required by modules 5.4 (RGI BWT) and other RGI-based AMR detection.

```shell
mkdir -p databases/RGI_databases
cd databases/RGI_databases

# Activate RGI environment
source rgi_env/bin/activate

# Download and load the CARD database
wget https://card.mcmaster.ca/latest/data
tar -xjf data ./card.json
rgi load --card_json card.json --local

# For RGI BWT, also load the CARD variants and WildCARD data
wget https://card.mcmaster.ca/latest/variants
tar -xjf variants
rgi card_annotation -i card.json
rgi load -i card.json --card_annotation card_database_*.fasta --local
rgi load --wildcard_annotation wildcard_database_*.fasta --wildcard_index *.txt --card_annotation card_database_*.fasta --local

deactivate
```

### ShortBRED Markers

Required by module 5.3 for marker-based AMR quantification.

```shell
mkdir -p databases/ShortBRED
# ShortBRED CARD markers can be generated from the CARD database using shortbred_identify.py
# or obtained from the CARD downloads page
```

### NCBI BLAST Database

Required by AMR detection modules for similarity searches.

```shell
mkdir -p databases/blast/db
cd databases/blast/db
update_blastdb.pl --decompress nt nr
```

### MetaPhlAn4 Database

Required by module 2.3 (optional).

```shell
conda activate metawrap-env  # or an environment with MetaPhlAn4
metaphlan --install --bowtie2db databases/metaphlan4_db
conda deactivate
```

---

## 3. Configuration

After installing all tools and databases, configure the pipeline:

1. Copy and edit the appropriate `.config` file in `pipelineScripts/configs/`:
   - Set `ROOT` to the full path of the `sequence_processing/` directory
   - Set `DATA_ROOT` and `PROCESSED_ROOT` to your data directories
   - Set `seqPath` to the full path of your raw sequence files
   - Set `sampleList` to the path of your sample ID manifest
   - Set `R1_ext` and `R2_ext` to match your raw sequence file naming convention
   - Update all database paths to match your installation locations
   - Update Conda environment paths to match your system

2. Update the metaWRAP config file (`which config-metawrap`) with:
   - `KRAKEN2_DB` — path to your Kraken2 database
   - `BMTAGGER_DB` — (only if using metaWRAP's built-in host removal; this pipeline uses Bowtie2 instead)

### Sample Manifest

Create a text file with one sample ID per line in `pipelineScripts/sample_lists/`.

The pipeline assumes sequence files are named:
`<sample_id>_<index_barcode>_S<sample_number>_L<lane_number>_R1_001.fastq.gz`

The sample IDs in the manifest should not include the paired-end extensions — only the common prefix up to the lane number (e.g., `<sample_id>_<index_barcode>_S<sample_number>_L<lane_number>`). Each sample needs only one entry.
