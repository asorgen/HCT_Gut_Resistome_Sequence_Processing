# HCT Gut Resistome — Sequence Processing Pipeline

This repository contains a modular metagenomic pipeline for processing HCT patient gut microbiome samples on the HPC (SLURM). It supports short-read, long-read, and hybrid assembly workflows.

---

## Repository Layout

```
HCT_Gut_Resistome_Sequence_Processing/
├── pipelineScripts/
│   ├── configs/                  # Config files per pipeline variant
│   │   ├── private.config        # Machine-specific paths (gitignored — you must create this)
│   │   ├── private.config.example
│   │   ├── functions.sh          # Core SLURM submission logic (run_module, module_setup, etc.)
│   │   ├── function_commands.sh  # Interactive helper functions (run_pipeline, check_jobs)
│   │   ├── Duke_short-read.config
│   │   ├── Duke_long-read.config
│   │   ├── Duke_hybrid-read.config
│   │   └── UNC_short-read.config
│   ├── pipeline_wrappers/        # Entry points — one per cohort/read-type
│   ├── sample_lists/             # One sample ID per line (header: #SampleID)
│   └── scripts/
│       ├── short_scripts/        # Modules 0.x–6 for Illumina reads
│       ├── long_scripts/         # Modules for Nanopore reads
│       ├── hybrid_scripts/       # Hybrid assembly modules
│       ├── helper_scripts/       # Utilities: parsing, DB download, gene processing
│       └── post_scripts/         # Summary and aggregation scripts
└── manuals/                      # Per-module documentation
```

---

## How the Pipeline Works

1. **Entry point**: `pipeline_wrappers/<dataset>_pipeline.sh`
   - Sources the corresponding `.config` file (which sources `private.config`)
   - Iterates through each sample in the `sample_lists/` file
   - Calls `run_module` for each enabled module, chaining SLURM dependencies automatically

2. **Module completion flags**: Each module writes a `COMPLETE/<ID>` flag file when it finishes. `run_module` checks for this flag before submitting — already-completed modules are skipped automatically.

3. **Module toggling**: Each module has a `run_<module>=true/false` flag in the `.config` file. Set to `false` to skip.

4. **Exclude samples**: Add sample IDs to `exclude_ids` in `private.config` to skip them in all wrappers.

---

## Configuration — private.config

`private.config` is gitignored and must be created on each machine from the example:

```bash
cp pipelineScripts/configs/private.config.example pipelineScripts/configs/private.config
```

Key variables to set:

| Variable | Description |
|---|---|
| `HPC_USER` | Your HPC username |
| `HPC_HOME` | Your home directory (e.g. `/users/asorgen`) |
| `HPC_PROJECTS` | Your lab project directory (e.g. `/projects/afodor_research3/asorgen`) |
| `HPC_SCRATCH` | Your scratch directory |
| `email_short` | Email for short-read pipeline SLURM notifications |
| `email_long` | Email for long/hybrid-read pipeline SLURM notifications |
| `SEQ_PATH_DUKE_SHORT` | Path to Duke Illumina raw FASTQ directory |
| `SEQ_PATH_UNC_SHORT` | Path to UNC Illumina raw FASTQ directory |
| `SEQ_PATH_DUKE_LONG` | Path to Duke Nanopore raw FASTQ directory |
| `exclude_ids` | Space-separated sample IDs to skip (e.g. failed QC samples) |

All config files and scripts derive their paths from `${HPC_PROJECTS}` — never hardcode personal paths.

---

## Sample Lists

Sample lists live in `pipelineScripts/sample_lists/`. Format:

```
#SampleID
SAMPLE001_GGACTCCT-ACTGCATA_S221_L004
SAMPLE002_CTCTCTAC-ACTGCATA_S134_L003
```

- The pipeline wrapper reads the full sequencing filename (used to locate FASTQ files)
- The sample ID (`ID`) is extracted as everything before the first `_`
- **Never hardcode arrays of sample IDs in scripts** — always read from the sample list file

---

## Running the Pipeline

```bash
# Run from the sequence_processing working directory on the HPC
cd /path/to/sequence_processing

# Short-read (Duke cohort)
nohup bash pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh >> Duke_short/LOGs/Duke_short_pipeline.out 2>&1 &

# Long-read
nohup bash pipelineScripts/pipeline_wrappers/Duke_long_pipeline.sh >> Duke_long/LOGs/Duke_long_pipeline.out 2>&1 &

# Stop after a specific module (pass the pipeline tag as arg)
bash pipelineScripts/pipeline_wrappers/Duke_short_pipeline.sh trim
```

### Restarting a module
Delete the relevant `COMPLETE/<ID>` flag file for any sample you want to rerun. The next pipeline run will re-submit that module for those samples.

### Checking progress
```bash
# Count completed samples for a module
squeue -u $USER
ls Duke_short/5.6_AA_amr_bins/COMPLETE/ | wc -l
```

---

## Post-Pipeline Summarization

After modules complete, aggregate results with `pipeline_summary.sh`:

```bash
# All modules
sbatch pipelineScripts/scripts/post_scripts/pipeline_summary.sh -p Duke_short

# Specific modules
bash pipelineScripts/scripts/post_scripts/pipeline_summary.sh -p Duke_short -m rgi_bwt,asm_gene_profiling
```

Output tables go to `<PROCESSED_ROOT>/Duke_short_tables/`.

---

## Rules for Modifying Scripts

- **Never hardcode personal paths** (`/projects/afodor_research3/asorgen/...`). Use `${HPC_PROJECTS}` sourced from `private.config`.
- **Never hardcode sample ID arrays** in scripts. Read IDs from the sample list file as the pipeline wrappers do.
- **Never hardcode email addresses** in `#SBATCH` headers. Use `${email_short}` or `${email_long}` from `private.config` in the opts strings, or omit `--mail-user` for standalone scripts.
- Scripts that need `${HPC_PROJECTS}` and aren't sourced through a config must source `private.config` themselves:
  ```bash
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  source "${SCRIPT_DIR}/../../configs/private.config"
  ```
- Python scripts with interactive fallback defaults should use generic placeholder paths (`/path/to/scratch/<SAMPLE_ID>...`), not personal paths.

---

## Active Work — HCT Study

This pipeline is being used to process HCT patient gut microbiome samples for the [[HCT Gut Resistome Study]]. The primary output feeding into downstream statistical analysis is from **Module 5.6** (`AA_amr_bins` — protein-level AMR detection from MAGs), producing per-sample gene count tables consumed by `HCT_Gut_Resistome_Heston_Analysis`.

Currently running: `D20248D1` 5.6_AA_amr_bins (in progress on HPC as of 2026-04-01).
