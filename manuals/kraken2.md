# Taxonomic classification with Kraken2

**Description**: Running Kraken2 on the reads will give us an idea of the taxonomic composition of the community. Running Kraken2 on the assembly will give us an idea what taxonomic groups were assembled better than others (the assembly process is heavily biased and should no be used to infer overall community composition).

## Module overview

  1. Runs Kraken2 on reads
  2. Runs Kraken2 on the assembly
  3. Runs Kraken-translate on the output
  4. Generates a kronagram of all files

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

## Run module

Assuming you have this github cloned, you will need to enter the root cloned directory and enter the following command:

### Short-read analysis

Run command
```shell
nohup sh ./SHORT/short_pipeline.sh kraken2 > SHORT/short_pipeline.out 2>&1 &
```

**Input:**
 - `evaluation/<SampleID>_final_assembly.fasta`
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`


### Long-read analysis

Run command
```shell
nohup sh ./LONG/long_pipeline.sh kraken2 > LONG/long_pipeline.out 2>&1 &
```

**Input:**
 - `evaluation/<SampleID>_final_assembly.fasta`
 - `clean_reads/<SampleID>_ont.fastq`


### Hybrid assembly analysis

Run command
```shell
nohup sh ./HYBRID/hybrid_pipeline.sh kraken2 > HYBRID/hybrid_pipeline.out 2>&1 &
```

**Input:**
 - `evaluation/<SampleID>_final_assembly.fasta`
 - `clean_reads/<SampleID>_1.fastq`
 - `clean_reads/<SampleID>_2.fastq`
 - `clean_reads/<SampleID>_ont.fastq`

## Output

- `kraken2/<SampleID>/*.krona` - files summarize taxonomy statistics to be fed into KronaTools
- `kraken2/<SampleID>/*.kraken2` - files contain the Kraken-estimated taxonomy of each read or contig
- `kraken2/<SampleID>/kronagram.html` - file contains all the taxonomy information from all the sample and assembly



