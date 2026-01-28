# Setup

## Required tools

 **metaWRAP**: MetaWRAP aims to be an **easy-to-use metagenomic wrapper suite** that accomplishes the core tasks of metagenomic analysis from start to finish: read quality control, assembly, visualization, taxonomic profiling, extracting draft genomes (binning), and functional annotation.

 **bbmap**: Download BBTools from [Sourceforge](https://sourceforge.net/projects/bbmap/)

## Required databases
 
 ### Indexed host genome
  **Host genome index for bmtagger**: you will need the bmtagger hg38 index to remove the human reads - see the metaWRAP database installation instructions. 

  To download and merge the human genome hg38
  ```shell
  mkdir BMTAGGER_INDEX
  cd BMTAGGER_INDEX
  wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/*fa.gz
  gunzip *fa.gz
  cat *fa > hg38.fa
  rm chr*.fa
  ```

  Install bmtool & srprism
  ```shell
  conda install bioconda::bmtool
  conda install bioconda::srprism
  ```

  Index the human genome. 
  *Note:* the file names of the indices must be exactly as specified for metaWRAP to recognize them.
  *Note:* indexing takes considerable memory and time (here - pass 100GB of RAM as -M parameter)
  ```shell
  bmtool -d hg38.fa -o hg38.bitmask
  srprism mkindex -i hg38.fa -o hg38.srprism -M 100000
  ```

  *Note:* metaWRAP looks for files hg38.bitmask and hg38.srprism
  **Don't forget to specify the BMTAGGER_DB variable in the config-metawrap file.** Run `which config-metawrap` to find it. 

 ### Kraken2 taxonomic database
  This database is custom-made, but any regular default Kraken2 database will do.

  **Don't forget to specify the KRAKEN2_DB variable in the config-metawrap file.** Run `which config-metawrap` to find it.

 ### CheckM database

## Short-read analysis
 To run this pipeline or any of the individual modules, edit `SHORT/short-read.config` such that 
 - `ROOT` indicates the full path to your cloned **HCT_Multi_Assembly_Approach_Analysis** directory
 - `seqPath` indicates the full path to your short-read illumina sequences
 - `sampleList` indicates the path to your sample ID manifest
 - `R1_ext` and `R2_ext` to reflect the precise sequence file naming convention of your raw sequence files.
 
 ### Sample manifest
 *This program runs under the assumption that sequence files are named 
 <sample_id>\_<index_barcode>\_S<sample_number>\_L<lane_number>\_R1_001.fastq.gz and <sample_id>\_<index_barcode>\_S<sample_number>\_L<lane_number>\_R1_001.fastq.gz.* 
 
 The sample IDs within the manifest do not need to include paired-end extensions, only <sample_id>\_<index_barcode>\_S<sample_number>\_L<lane_number>. And thus only require once instance of the sample ID per sample.

