#! /usr/bin/env python3

# Usage
# python3 amr_count.py -f AMR_directory -s LONG

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Script to count AMR in files.")
parser.add_argument(
    "--folder", "-f",
    type=str,
    required=True,
    help="Path to the directory containing files to process"
)
parser.add_argument(
    "--suffix", "-s",
    type=str,
    required=True,
    help="Suffix or label to append to output file."
)

args = parser.parse_args()

directory = args.folder
label = args.suffix

contigs = [f for f in os.listdir(directory) if f.endswith('rgi.txt')]

dc = {}
for name in contigs:
    path = directory + "/" + name
    sample = name.replace(".rgi.txt", "")
    # print(f"Path: {path}")
    dc[sample] = pd.read_csv(path, sep="\t")
    # print(dc)

# Initialize an empty DataFrame
dfc = pd.DataFrame()

# Step 1: Collect all unique 'Best_Hit_ARO' values across all files
unique_gene_symbols = set()
for name in dc:
    unique_gene_symbols.update(dc[name]['Best_Hit_ARO'].unique())

# Convert the set to a sorted list (to ensure order)
unique_gene_symbols = sorted(list(unique_gene_symbols))

# Step 2: Create a DataFrame with all unique 'Best_Hit_ARO' as the index
dfc = pd.DataFrame(index=unique_gene_symbols)

# Step 3: Populate 'dfc' with counts from each file
for name in dc:
    counts = dc[name]['Best_Hit_ARO'].value_counts(sort=False)
    dfc[name] = counts  # Aligns based on the 'Best_Hit_ARO' index

dfc = dfc.fillna(0).astype(int)
# print(dfc)
outputFile = directory + "/" + "RGI_" + label + ".tsv"
dfc.to_csv(outputFile, sep = '\t')
