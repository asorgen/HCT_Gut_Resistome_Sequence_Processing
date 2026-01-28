#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 12:58:28 2025

@author: aliciasorgen
"""

import pandas as pd
import argparse
import os

def read_stats_file(file_path):
    """
    Reads a .stats file into a pandas DataFrame.
    Assumes a tab-delimited file format. Modify `sep` if needed.
    """
    try:
        column_names=['bin', 'bin_classification']
        df = pd.read_csv(file_path, sep="\t", engine='python', header=None, names=column_names)  # Change separator if necessary
        print("File loaded successfully!")
        return df
    except Exception as e:
        return None

def append_sample_id(df, sample_id):
    """
    Adds a 'SampleID' column to the DataFrame with all entries set to the provided sample_id string.
    """
    df["SampleID"] = sample_id
    return df

def split_bin_column(df):
    """
    Splits the 'bin' column at the second occurrence of '.' and creates two new columns:
    - 'bin' contains 'bin.*'
    - 'bin_type' contains everything after the second '.'
    """
    if "bin" not in df.columns:
        print("Column 'bin' not found in the dataset.")
        return df

    # Splitting the 'bin' column at the second '.'
    split_bins = df["bin"].str.split(".", n=2, expand=True)

    # Assign the split values to new columns
    df["bin"] = split_bins[0] + "." + split_bins[1]  # First two parts
    df["bin_type"] = split_bins[2]  # Remaining part
    df["bin_type"] = df["bin_type"].str.replace(".fa", "", regex=False)

    return df

def process_and_append(file_path, sample_id, master_file):
    """
    Processes a .stats file, updates it, and appends it to the master file.
    """
    df = read_stats_file(file_path)
    
    if df is not None:
        # Process the DataFrame
        df = append_sample_id(df, sample_id)
        df = split_bin_column(df)

        # Append to master file
        write_header = not os.path.exists(master_file)  # Only write header if file doesn't exist
        df.to_csv(master_file, sep="\t", index=False, mode="a", header=write_header)

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

