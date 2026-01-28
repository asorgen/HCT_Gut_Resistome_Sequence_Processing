# Assembly refinement and evaluation


**Description**:

## Module overview
 1. filter contigs to a min length of 1.5kb
 2. evaluate assembly statistics 


## Setup

See [setup](manuals/setup.md) for requirements. 


## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

**Input:**
 - `assembly/<SampleID>/assembly.fasta`

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh evaluation > SHORT/short_pipeline.out 2>&1 &
```

### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh evaluation > LONG/long_pipeline.out 2>&1 &
```

### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh evaluation > HYBRID/hybrid_pipeline.out 2>&1 &
```


## Output Overview

- `evaluation/<SampleID>_draft_assembly.fasta` - original assembly taken straight from the assembly module
- `evaluation/<SampleID>_final_assembly.fasta` - filtered assembly containing reads >1.5 kb
- `evaluation/<read_type>_read_draft_assembly_stats.txt` - summary of all draft assemblies from pipeline
- `evaluation/<read_type>_read_final_assembly_stats.txt` - summary of all final assemblies from pipeline

*Note*: the assembly stats will not be generated until all samples from the pipeline have been processed
