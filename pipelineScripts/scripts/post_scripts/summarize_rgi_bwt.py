#! /usr/bin/env python3

# Usage
# python3 summarize_rgi_bwt.py -f $bwt_dir -l $pipeline


# %%
import os
import pandas as pd
import argparse
# import numpy as np


# %% Argument parser
parser = argparse.ArgumentParser(description="Script to count AMR in files.")
parser.add_argument(
    "--folder", "-f",
    type=str,
    required=True,
    help="Path to the directory containing files to process"
)
parser.add_argument(
    "--label", "-l",
    type=str,
    required=True,
    help="Label to append to output file."
)
parser.add_argument(
    "--reads", "-r",
    type=str,
    required=True,
    help="File containing read summary."
)
parser.add_argument(
    "--mapped", "-m",
    type=int,
    required=True,
    help="Minimum number of mapped reads."
)
parser.add_argument(
    "--mapq", "-M",
    type=int,
    required=True,
    help="Minimum MAPQ score."
)
parser.add_argument(
    "--coverage", "-c",
    type=int,
    required=True,
    help="Minimum percent coverage."
)
parser.add_argument(
    "--outdir", "-o",
    type=str,
    required=False,
    default=None,
    help="Output directory. If not provided, derived from folder path."
)

# %% Naming variables
# Try to parse args, set defaults if running interactively
try:
    args = parser.parse_args()
    # Name variables from command-line arguments
    directory = args.folder
    pipeline = args.label
    fastqSum = args.reads
    min_reads = args.mapped
    min_mapq = args.mapq
    min_cov = args.coverage
    outdir = args.outdir
except SystemExit:
    # Set default values for interactive development
    print("Running in interactive mode - using default values")
    directory = "/Users/aliciasorgen/scratch/kma_output"
    pipeline = "Duke_short"
    fastqSum = "/Users/aliciasorgen/scratch/Duke_short_post-QC_report.tsv"
    min_reads = 5
    min_mapq = 0
    min_cov = 0
    outdir = None

alignment_folder = directory.split("/")[-1]
aligner = alignment_folder.split("_")[0]
if outdir is None:
    outdir = os.path.dirname(os.path.dirname(directory)) + "/" + pipeline + "_tables"

# %% Importing and summarizing fastq files

fastq = pd.read_csv(fastqSum, sep="\t")
fastq = fastq[['Sample', 'FastQC_mqc-generalstats-fastqc-avg_sequence_length', 'FastQC_mqc-generalstats-fastqc-total_sequences']]
fastq = fastq.rename(columns={'FastQC_mqc-generalstats-fastqc-avg_sequence_length': 'avg_sequence_length',
                              'FastQC_mqc-generalstats-fastqc-total_sequences': 'total_sequences'})
fastq['Sample'] = fastq['Sample'].str[:-2]

summary_df = fastq.groupby('Sample').agg({
    'avg_sequence_length': 'mean',
    'total_sequences': 'sum'
}).reset_index()

summary_df['avg_sequence_length'] = summary_df['avg_sequence_length'].round(0).astype(int)

# %% Setting file names for importing

file_ext = ".rgi_" + aligner + ".txt"

# Save all the RGI bwt output files within a list ('files')
files = [f for f in os.listdir(directory) if f.endswith(file_ext)]


# %% Looping through files and calculating FPKM

dc = {}
df_list = []
hits_df = pd.DataFrame(columns=['SampleID', 'Num_AMR_Hits'])
for name in files:
    path = directory + "/" + name
    sample = name.replace(file_ext, "")
    df = pd.read_csv(path, sep="\t")
    df['SampleID'] = sample
    index = summary_df[summary_df['Sample'] == sample].index
    R = summary_df.loc[index, 'avg_sequence_length'].values[0]
    N = summary_df.loc[index, 'total_sequences'].values[0] /2
    
    # Filter amrfinder to only keep 'Average Percent Coverage' >= 50%
    df = df[(df['Average Percent Coverage'] >= min_cov) &
            (df['Completely Mapped Reads'] > min_reads) &
            (df['Average MAPQ (Completely Mapped Reads)'] >= min_mapq)]
    
    df['ARO Accession'] = 'ARO_' + df['ARO Accession'].astype(str)
    df['Reference Length'] = pd.to_numeric(df['Reference Length'], errors='coerce')
    
    # Calculate FPKM    
    df['FPKM'] = (df['All Mapped Reads'] * 10**9) / (df['Reference Length'] * N)
    
    dc[sample] = df
    df_list.append(df)
    
    # Append to hits_df
    new_row = {'SampleID': sample, 'Num_AMR_Hits': len(df['ARO Accession'])}
    hits_df = pd.concat([hits_df, pd.DataFrame([new_row])], ignore_index=True)
    

# %% Concatenate all at once (much faster than repeated appending)
full_df = pd.concat(df_list, ignore_index=True)

# %% Create abundance table
# Initialize an empty DataFrame
dfc = pd.DataFrame()

# Collect all unique 'ARO Accession' values across all files
unique_ARO = set()
for name in dc:
    unique_ARO.update(dc[name]['ARO Accession'].unique())

# Convert the set to a sorted list (to ensure order)
unique_ARO = sorted(list(unique_ARO))

# Create a DataFrame with all unique 'ARO Accession' as a column
dfc = pd.DataFrame({'ARO Accession': unique_ARO})

# Populate 'dfc' with counts from each file
for name in dc:
    counts = dc[name].set_index('ARO Accession').reindex(dfc['ARO Accession'])['FPKM'].values
    dfc[name] = counts  # Aligns based on the 'ARO Accession' index

# dfc = dfc.fillna(0).astype(float)

# fill NAs with 0 and convert only columns from 2nd onward
dfc.iloc[:, 1:] = dfc.iloc[:, 1:].fillna(0).astype(float)

# %% Set output names
gene_table_name = "_rgi_bwt_" + aligner + "_FPKM_R" + str(min_reads) + "_M" + str(min_mapq) + "_C" + str(min_cov) + ".tsv"
concat_output_name = "_rgi_bwt_" + aligner + "_output_R" + str(min_reads) + "_M" + str(min_mapq) + "_C" + str(min_cov) + ".tsv"
hits_name = "_rgi_bwt_" + aligner + "_hits_R" + str(min_reads) + "_M" + str(min_mapq) + "_C" + str(min_cov) + ".tsv"

table_outputFile = outdir + "/" + pipeline + gene_table_name
concat_outputFile = outdir + "/" + pipeline + concat_output_name
hits_outputFile = outdir + "/" + pipeline + hits_name



# %% Save outputs
dfc.to_csv(table_outputFile, sep = '\t')
full_df.to_csv(concat_outputFile, sep = '\t')
hits_df.to_csv(hits_outputFile, sep = '\t')
print(hits_outputFile)
