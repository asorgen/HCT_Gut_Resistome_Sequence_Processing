# Functional annotation of bins

**Description:** <add something>

## Module overview

<add something>

## Setup

See [setup](manuals/setup.md) for requirements. 


## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

**Input:**
- `reassemble_bins/reassembled_bins/`

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh annotation > SHORT/short_pipeline.out 2>&1 &
```


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh annotation > LONG/long_pipeline.out 2>&1 &
```


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh annotation > HYBRID/hybrid_pipeline.out 2>&1 &
```


## Output Overview

- `annotate_bins/bin_funct_annotations/` - contains the functional annotations of each bin in .GFF format