#! /usr/bin/env python3

# Usage
# python3 summarize_shortbred.py -f $sb_dir -l $pipeline

import os
import pandas as pd
import argparse
import numpy as np

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
    help="File containing read summary data."
)

args = parser.parse_args()

directory = args.folder
pipeline = args.label
fastqSum = args.reads

if not os.path.isdir(directory):
    print(f"ShortBRED directory not found: {directory}. Skipping...")
    import sys; sys.exit(0)

# %%

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

# %%

file_ext = '.shortbred.tsv'

files = [f for f in os.listdir(directory) if f.endswith(file_ext)]
# print(files)

dc = {}
for name in files:
    path = directory + "/" + name
    sample = name.replace(file_ext, "")
    df = pd.read_csv(path, sep="\t")
    df[['X', 'GenBank_Acc', 'ARO_Acc', 'Annotation']] = df['Family'].str.split('|', expand=True)
    index = summary_df[summary_df['Sample'] == sample].index
    # print(index)
    R = summary_df.loc[index, 'avg_sequence_length'].values[0]
    N = summary_df.loc[index, 'total_sequences'].values[0]
    df['L'] = df['TotMarkerLength'] * 3
    df['L_prime'] = np.where(df['L'] > (0.95*R),
                             df['L'] - (0.9 * R) + 1,   # if
                             R - df['L'] - 1)          # else    
    df['RPKM'] = df['Hits'] / ((df['L_prime'] / 1000) * (N / 1000000))
    dc[sample] = df


# %%
# Initialize an empty DataFrame
dfc = pd.DataFrame()

# Collect all unique 'ARO_Acc' values across all files
unique_ARO = set()
for name in dc:
    unique_ARO.update(dc[name]['ARO_Acc'].unique())

# Convert the set to a sorted list (to ensure order)
unique_ARO = sorted(list(unique_ARO))

# Create a DataFrame with all unique 'ARO_Acc' as a column
dfc = pd.DataFrame({'ARO_Acc': unique_ARO})

# Populate 'dfc' with counts from each file
for name in dc:
    counts = dc[name].set_index('ARO_Acc').reindex(dfc['ARO_Acc'])['RPKM'].values
    dfc[name] = counts  # Aligns based on the 'ARO_Acc' index

outputFile = directory + "/" + pipeline +  "_shortbred.tsv"
dfc.to_csv(outputFile, sep = '\t')