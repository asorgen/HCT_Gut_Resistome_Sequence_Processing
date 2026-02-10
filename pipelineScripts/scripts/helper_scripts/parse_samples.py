#! /usr/bin/env python3

# Usage
# python3 parse_samples.py -c Duke -f Samples.txt -o parsed_samples.tsv
# %% Import libraries
import re
import pandas as pd
import argparse

# %% Parse arguments
parser = argparse.ArgumentParser(description="Parse sample IDs from Samples.txt into structured components.")
parser.add_argument(
    "--cohort", "-c",
    type=str,
    required=True,
    help="Sample cohort"
)
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

# %% Set pattern
if args.cohort == "Duke":
    pattern = re.compile(r'^D(\d+)(D-?\d+|PRE|NPE)_([ACGT]+)-([ACGT]+)_S(\d+)_L(\d+)$')
elif args.cohort == "UNC":
    pattern1 = re.compile(r'^BMT(\d+)(D-?\d+|d-?\d+|PRE|NPE|pre)_([ACGT]+)-([ACGT]+)_S(\d+)_L(\d+)$')
    pattern2 = re.compile(r'^BMT(\d+)(D-?\d+|d-?\d+|PRE|NPE|pre)_([ACGT]+)-([ACGT]+)_S(\d+)$')

records = []
with open(args.file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if args.cohort == "Duke":
            m = pattern.match(line)
        elif args.cohort == "UNC":
            m = pattern1.match(line)
            m2 = pattern2.match(line)
        if m:
            if args.cohort == "Duke":
                patient_id = "D" + m.group(1)
            elif args.cohort == "UNC":
                patient_id = "BMT" + m.group(1)
            timepoint_raw = m.group(2)
            sample_id = patient_id + timepoint_raw
            # Strip leading 'D' from numeric timepoints
            if timepoint_raw.startswith('D'):
                timepoint = timepoint_raw[1:]
            elif timepoint_raw.startswith('d'):
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
        elif m2:
            patient_id = "BMT" + m2.group(1)
            timepoint_raw = m2.group(2)
            sample_id = patient_id + timepoint_raw
            # Strip leading 'D' from numeric timepoints
            if timepoint_raw.startswith('D'):
                timepoint = timepoint_raw[1:]
            elif timepoint_raw.startswith('d'):
                timepoint = timepoint_raw[1:]
            else:
                timepoint = timepoint_raw
            records.append({
                'IlluminaID': line,
                'SampleID': sample_id,
                'PatientID': patient_id,
                'Timepoint': timepoint,
                'Barcode1': m2.group(3),
                'Barcode2': m2.group(4),
                'SampleNumber': int(m2.group(5)),
                'Lane': '',
            })
        else:
            print(f"WARNING: Could not parse line: {line}")

df = pd.DataFrame(records)

if args.output:
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote {len(df)} samples to {args.output}")
else:
    print(df.to_string(index=False))
