# Taxonomic bin classification

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
nohup sh ./SHORT/short_pipeline.sh classify > SHORT/short_pipeline.out 2>&1 &
```


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh classify > LONG/long_pipeline.out 2>&1 &
```


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh classify > HYBRID/hybrid_pipeline.out 2>&1 &
```


## Output Overview

- `classify_bins/bin_taxonomy.tab` - taxonomic annotation