#!/usr/bin/env python2.7

import argparse

def filter_fastq(input_file, output_file):
    seen_headers = set()  # Set to store seen headers

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        while True:
            try:
                header = infile.readline().strip()
                sequence = infile.readline().strip()
                separator = infile.readline().strip()
                quality = infile.readline().strip()

                # Check if the end of the file is reached
                if not header or not sequence or not separator or not quality:
                    break

                # Skip the read if the header has already been encountered
                if header in seen_headers:
                    continue

                # Add the header to the set
                seen_headers.add(header)

                # Skip the read if the sequence contains invalid characters
                if not all(char in "ATCGN" for char in sequence):
                    continue

                # Write valid reads to the output file
                outfile.write("{0}\n{1}\n{2}\n{3}\n".format(header, sequence, separator, quality))
            except StopIteration:
                break

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Filter reads from a FASTQ file based on sequence validity and unique headers.")
    parser.add_argument("input", help="Input FASTQ file")
    parser.add_argument("output", help="Output FASTQ file for filtered reads")

    # Parse the arguments
    args = parser.parse_args()

    # Call the filter function with the provided arguments
    filter_fastq(args.input, args.output)

if __name__ == "__main__":
    main()