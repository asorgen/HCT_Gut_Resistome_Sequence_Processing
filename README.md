# HCT Multi-Assembly Approach Analysis
Computational workflows for metagenomic tasks for short-reads, long-reads, and hybrid assemblies packaged with metaWRAP and Conda.

### Table of Contents

 1. [Setup](manuals/setup.md)
 2. [Running a module](manuals/run_module.md)
 3. Available modules:
    - [**Pre-processing** metagenomic data](manuals/preprocessing.md)
    - [Metagenomic **assembly**](manuals/assembly.md)
    - [Assembly refinement and **evaluation**](manuals/evaluation.md)
    - [Taxonomic classification with **Kraken2**](manuals/kraken2.md)
    - [Metagenomic **binning**](manuals/binning.md)
    - [Metagenomic bin **refinement**](manuals/refine.md)
    - [Metagenomic bin **reassembly**](manuals/reassembly.md)
    - [Taxonomic bin **classification**](manuals/classify_bins.md)
    - [**Functional annotation** of bins](manuals/annotation.md)
    - [**AMR detection**](manuals/amr_search.md)




### Quickstart
If you are working on a Slurm workload manager, this command is an example of how to run the complete pipeline. Otherwise, you will need to change options.
```shell
nohup sh ./SHORT/scripts/short_pipeline.sh > SHORT/short_pipeline.out 2>&1 &
nohup sh ./LONG/scripts/long_pipeline.sh > LONG/long_pipeline.out 2>&1 &
nohup sh ./HYBRID/scripts/hybrid_pipeline.sh > HYBRID/hybrid_pipeline.out 2>&1 &
```