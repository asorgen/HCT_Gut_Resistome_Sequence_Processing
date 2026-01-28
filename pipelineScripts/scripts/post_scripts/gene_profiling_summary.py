#!/usr/bin/env python3
import pandas as pd
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a .stats file and append a SampleID column.")
    parser.add_argument("file_path", type=str, help="Path to the .stats file")
    parser.add_argument("sample_id", type=str, help="String to fill the SampleID column")
    parser.add_argument("master_file", type=str, help="Path to the master output file")
    
    args = parser.parse_args()
    
    # Read in arguments
    file_path = args.file_path
    sample_id = args.sample_id
    master_file = args.master_file
    
    # Read in input file
    df = pd.read_csv(file_path, sep="\t", engine='python')  # Change separator if necessary
    
    # Add a 'SampleID' column to the DataFrame with all entries set to the provided sample_id string.
    df["SampleID"] = sample_id
    
    # Remove Bakta annotations
    # df = df[df['annotated_with'] != 'Bakta']

    # Append to master file
    write_header = not os.path.exists(master_file)  # Only write header if file doesn't exist
    df.to_csv(master_file, sep="\t", index=False, mode="a", header=write_header)