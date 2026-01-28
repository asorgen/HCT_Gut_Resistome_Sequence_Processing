#! /usr/bin/env python3

import pandas as pd
import argparse
import os
import numpy as np

def read_stats_file(file_path):
    """
    Reads a .stats file into a pandas DataFrame.
    Assumes a tab-delimited file format. Modify `sep` if needed.
    """
    try:
        column_names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        df = pd.read_csv(file_path, sep="\t", engine='python', header=None, names=column_names)  # Change separator if necessary
        
        new_path = file_path.split("/")
        file = new_path[-1]
        binner_name = file.replace(".gff", "")
        split_binner_name = binner_name.split(".")
        df["bin"] = split_binner_name[0] + "." + split_binner_name[1]
        df["bin_type"] = split_binner_name[2]
        
        return df
    except Exception as e:
        # print(f"Error loading file: {e}")
        return None

def append_sample_id(df, sample_id):
    """
    Adds a 'SampleID' column to the DataFrame with all entries set to the provided sample_id string.
    """
    df["SampleID"] = sample_id

    return df

def extract_by_prefix(row, prefix):
    for cell in row:
        if isinstance(cell, str) and cell.startswith(prefix):
            return cell
    return np.nan

def split_attributes_column(df):
    split_attributes = df["attributes"].str.split(";", expand=True)

    df['ID'] = split_attributes.apply(extract_by_prefix, axis=1, prefix='ID=')
    df['gene'] = split_attributes.apply(extract_by_prefix, axis=1, prefix='gene=')
    df['inference'] = split_attributes.apply(extract_by_prefix, axis=1, prefix='inference=')
    df['locus_tag'] = split_attributes.apply(extract_by_prefix, axis=1, prefix='locus_tag=')
    df['product'] = split_attributes.apply(extract_by_prefix, axis=1, prefix='product=')
    df['eC_number'] = split_attributes.apply(extract_by_prefix, axis=1, prefix='eC_number=')

    return df

def process_and_append(file_path, sample_id, master_file):
    df = read_stats_file(file_path)
    
    if df is not None:
        # Process the DataFrame
        df = append_sample_id(df, sample_id)
        df = split_attributes_column(df)

        # Append to master file
        write_header = not os.path.exists(master_file)  # Only write header if file doesn't exist
        df.to_csv(master_file, sep="\t", index=False, mode="a", header=write_header)
        # print(f"Appended processed data from {file_path} to {master_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a .stats file and append a SampleID column.")
    parser.add_argument("file_path", type=str, help="Path to the .stats file")
    parser.add_argument("sample_id", type=str, help="String to fill the SampleID column")
    parser.add_argument("master_file", type=str, help="Path to the master output file")
    
    args = parser.parse_args()
    
    file_path = args.file_path
    sample_id = args.sample_id
    master_file = args.master_file
    
    # Process and append data to the master file
    process_and_append(file_path, sample_id, master_file)

