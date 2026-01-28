#! /usr/bin/env python3

# Usage
# python3 amr_count.py -f AMR_directory -s LONG

import os
import pandas as pd
import argparse

# def main(folder, suffix):
#     # Example: Logic to handle files in the specified folder
#     print(f"Processing files in folder: {folder}")
#     # print(f"Using suffix: {suffix}")

#     # Example logic: List all files in the folder
#     if not os.path.isdir(folder):
#         print(f"The specified folder '{folder}' does not exist.")
#         return

#     for file_name in os.listdir(folder):
#         print(f"Processing file: {file_name}")

# if __name__ == "__main__":

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
# main(folder=args.folder, suffix=args.suffix)

directory = args.folder
label = args.suffix
file_ext = '.amrfinder.txt'

contigs = [f for f in os.listdir(directory) if f.endswith(file_ext)]
# print(contigs)

dc = {}
for name in contigs:
    path = directory + "/" + name
    sample = name.replace(file_ext, "")
    # print(f"Path: {path}")
    dc[sample] = pd.read_csv(path, sep="\t")
    # print(dc)

# Initialize an empty DataFrame
dfc = pd.DataFrame()

# Step 1: Collect all unique 'Gene symbol' values across all files
unique_gene_symbols = set()
# unique_accession = set()
for name in dc:
    unique_gene_symbols.update(dc[name]['Gene symbol'].unique())
    # unique_accession.update(dc[name]['Accession of closest sequence'].unique())

# # Convert the set to a sorted list (to ensure order)
# unique_accession = sorted(list(unique_accession))
unique_gene_symbols = sorted(list(unique_gene_symbols))

# Step 2: Create a DataFrame with all unique 'Gene symbol' as the index
# dfc = pd.DataFrame({'Accession': unique_accession})
# dfc = pd.DataFrame({'Annotation': unique_gene_symbols})
dfc = pd.DataFrame(index=unique_gene_symbols)

# Step 3: Populate 'dfc' with counts from each file
for name in dc:
    counts = dc[name]['Gene symbol'].value_counts(sort=False)
    dfc[name] = counts  # Aligns based on the 'Gene symbol' index

dfc = dfc.fillna(0).astype(int)
# print(dfc)
outputFile = directory + "/" + "AMRFinder_" + label + ".tsv"
dfc.to_csv(outputFile, sep = '\t')
