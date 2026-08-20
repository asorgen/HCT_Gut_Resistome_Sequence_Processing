#!/usr/bin/env python3

"""
Generates list of samples from paired-end reads.

Longer description if needed.
Usage: ./build_sample_list.py <directory> <output file>
"""

import re
import sys
from pathlib import Path

KNOWN_EXTENSIONS = {".gz", ".fastq", ".fq", ".fasta", ".fa"}
ILLUMINA_PAIRED_END_PATTERN = r"_R[12].*$"

def curate_file_list(dir_path):
    """Return a list of file Paths in dir_path, excluding subdirectories."""
    return [f for f in dir_path.iterdir() if f.is_file()]

def strip_file_extensions(list_of_files):
    """Return a list of file Paths stripped of file extensions."""
    new_list = []
    for f in list_of_files:
        while f.suffix in KNOWN_EXTENSIONS:
            f = f.with_suffix("")
        new_list.append(f)
    return new_list

def strip_sample_suffix(list_of_files):
    """Return a list of names with _R1/_R2 and everything after removed."""
    new_list = []
    for f in list_of_files:
        name = f.name
        cleaned = re.sub(ILLUMINA_PAIRED_END_PATTERN, "", name)
        new_list.append(cleaned)
    return new_list

def deduplicate(names):
    """Return a sorted list of unique names."""
    return sorted(set(names))

def write_sample_list(names, output_path):
    """Write names to output_path, one per line, with #SampleID header."""
    with open(output_path, "w") as fh:
        fh.write("#SampleID\n")
        for name in names:
            fh.write(f"{name}\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit("Usage: ./build_sample_list.py <directory> <output>")
    
    input_dir = Path(sys.argv[1])
    if not input_dir.is_dir():
        sys.exit(f"Not a valid directory: {input_dir}")

    output_file = Path(sys.argv[2])
    if output_file.exists():
        sys.exit(f"Output file already exists: {output_file}")        

    files = curate_file_list(input_dir)
    stripped_files = strip_file_extensions(files)
    sample_names = strip_sample_suffix(stripped_files)
    unique_sample_names = deduplicate(sample_names)
    
    print(len(unique_sample_names), "total samples")
    # print(unique_sample_names[0])
    write_sample_list(unique_sample_names, output_file)