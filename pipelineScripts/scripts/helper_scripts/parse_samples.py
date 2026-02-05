#! /usr/bin/env python3

# Usage
# python3 parse_samples.py -f Samples.txt -o parsed_samples.tsv
# %%
import re
import pandas as pd
import argparse

# %%
parser = argparse.ArgumentParser(description="Parse sample IDs from Samples.txt into structured components.")
parser.add_argument(
    "--file", "-f",
    type=str,
    required=True,
    help="Path to the Samples.txt file"
)
parser.add_argument(
    "--output", "-o",
    type=str,
    default=None,
    help="Path to output TSV file (optional; prints to stdout if not provided)"
)

args = parser.parse_args()

pattern = re.compile(r'^D(\d+)(D-?\d+|PRE|NPE)_([ACGT]+)-([ACGT]+)_S(\d+)_L(\d+)$')

records = []
with open(args.file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = pattern.match(line)
        if m:
            patient_id = "D" + m.group(1)
            timepoint_raw = m.group(2)
            sample_id = patient_id + timepoint_raw
            # Strip leading 'D' from numeric timepoints
            if timepoint_raw.startswith('D'):
                timepoint = timepoint_raw[1:]
            else:
                timepoint = timepoint_raw
            records.append({
                'IlluminaID': line,
                'SampleID': sample_id,
                'PatientID': patient_id,
                'Timepoint': timepoint,
                'Barcode1': m.group(3),
                'Barcode2': m.group(4),
                'SampleNumber': int(m.group(5)),
                'Lane': m.group(6),
            })
        else:
            print(f"WARNING: Could not parse line: {line}")

df = pd.DataFrame(records)

if args.output:
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(df)} samples to {args.output}")
else:
    print(df.to_string(index=False))
