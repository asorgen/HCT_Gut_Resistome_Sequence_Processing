# Metagenomic bin reassembly

**Description:**

## Module overview


## Setup

See [setup](manuals/setup.md) for requirements. 


## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh reassemble > SHORT/short_pipeline.out 2>&1 &
```

**Input:**
 - `refine_bins/<SampleID>/metawrap_70_10_bins/`
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh reassemble > LONG/long_pipeline.out 2>&1 &
```

**Input:**
 - `refine_bins/<SampleID>/metawrap_70_10_bins/`
 - `clean_reads/<SampleID>_ont.fastq`


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh reassemble > HYBRID/hybrid_pipeline.out 2>&1 &
```

**Input:**
 - `refine_bins/<SampleID>/metawrap_70_10_bins/`
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`
 - `clean_reads/<SampleID>_ont.fastq`



## Output Overview

- `reassemble_bins/reassembled_bins` - best quality bins after reassembly
- `reassemble_bins/original_bins.stats` - statistics for the original consolidated bins
- `reassemble_bins/reassembled_bins.png` - CheckM plot of the final bins
- `reassemble_bins/reassembled_bins.stats` - statistics for the best bin outputs
- `reassemble_bins/reassembly_results.png` - plot comparing the original and reassembled bin sets
