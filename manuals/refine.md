# Metagenomic bin refinement

**Description:** During refinement, metaWRAP will have to chose the best version of each bin between 7 different versions of each bin. It will dynamically adjust to prioritize the bin quality that you desire.

## Module overview

 1. Refine bins
 2. Run CheckM on refined bins
 3. Consolidate the refined bins
 4. Dereplicate the consolidated refined bins
 5. Run CheckM on the final (dereplicated, consolidated, refined) bins
 6. Plot completion and contamination rankings of final bins

## Setup

See [setup](manuals/setup.md) for requirements. 

## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

**Input:**
- `binning/concoct_bins`
- `binning/maxbin2_bins`
- `binning/metabat2_bins`


### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh refine > SHORT/short_pipeline.out 2>&1 &
```


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh refine > LONG/long_pipeline.out 2>&1 &
```


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh refine > HYBRID/hybrid_pipeline.out 2>&1 &
```


## Output Overview

- `refine_bins/<SampleID>/metawrap_70_10_bins/` - consolidated and refined bins for downstream analysis
- `refine_bins/<SampleID>/*.contigs` - contigs from various binning strategies
- `refine_bins/<SampleID>/*.stats` - statistics from each binning strategy


