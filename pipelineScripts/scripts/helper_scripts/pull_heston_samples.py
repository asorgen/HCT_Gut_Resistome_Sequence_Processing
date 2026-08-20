#!/usr/bin/env python3

"""
Pulls out the samples used in the Heston analysis from the entire SRA fastq-dump file list.

Longer description if needed.
Usage: ./pull_heston_samples.py <sample list> <metadata file> <output file>
"""

import re
import pandas as pd
import sys
from pathlib import Path
from build_sample_list import write_sample_list

SRA_PATTERN = r"_SRR.*$"

def strip_sra_id_suffix(list_of_files):
    """Return a list of names with the '_SRR*' removed."""
    new_list = []
    for f in list_of_files:
        name = f
        cleaned = re.sub(SRA_PATTERN, "", name)
        new_list.append(cleaned)
    return new_list

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit("Usage: ./pull_heston_samples.py <sample list> <metadata file> <output file>")
    
    sample_list = Path(sys.argv[1])
    if not sample_list.exists():
        sys.exit(f"This file could not be found: {sample_list}")
    
    metadata = Path(sys.argv[2])
    if not metadata.exists():
        sys.exit(f"This file could not be found: {metadata}")

    output_file = Path(sys.argv[3])
    if output_file.exists():
        sys.exit(f"Output file already exists: {output_file}")

    sample_df = pd.read_csv(sample_list)
    all_seq_ids = sample_df["#SampleID"]

    all_sample_ids = strip_sra_id_suffix(all_seq_ids)

    mapping = dict(zip(all_sample_ids, all_seq_ids))
    assert len(all_sample_ids) == len(all_seq_ids), \
    f"Length mismatch: {len(all_sample_ids)} vs {len(all_seq_ids)}"

    meta_df = pd.read_csv(metadata)
    keep_sample_ids = meta_df["SampleID"]

    kept_seq_ids = [mapping[sid] for sid in keep_sample_ids if sid in mapping]
    write_sample_list(kept_seq_ids, output_file)
