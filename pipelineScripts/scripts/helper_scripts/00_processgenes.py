#!/usr/bin/env python3
# %% Description -----------------------------------------------------------------------------------
#This script matches gene annotations in the following order of priority: AMRFinder, RGI, Bakta
'''
Slurm Usage Example: 
python scripts/helper_scripts/00_processgenes.py \
    $gene_list \ # asm_gene_annotation/D13004D15/D13004D15_geneslist.tsv
    $amr_output \ # asm_gene_annotation/D13004D15/D13004D15.amrfinder.tsv
    $rgi_output \ # asm_gene_annotation/D13004D15/D13004D15.rgi.txt
    $bakta_output \ # asm_gene_annotation/D13004D15/bakta/D13004D15.tsv
    $outputFile # asm_gene_annotation/D13004D15/D13004D15_gene_annotations.tsv
'''

'''
Test Usage: 
python3 scripts/helper_scripts/00_processgenes.py 
    /Users/aliciasorgen/Downloads/D13004D15_geneslist.tsv 
    /Users/aliciasorgen/Downloads/D13004D15.amrfinder.tsv 
    /Users/aliciasorgen/Downloads/D13004D15.rgi.txt
    /Users/aliciasorgen/Downloads/bakta/D13004D15.tsv 
    /Users/aliciasorgen/Downloads/D13004D15_gene_annotations.tsv
'''
# %% Import modules --------------------------------------------------------------------------------
import pandas as pd
import sys
import numpy as np
import argparse

# %% Argument parser -------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="This script merges and consolidates the AMR gene " \
"outputs from RGI, AMRFinder, and Bakta.")

parser.add_argument(
    "--genes",
    type=str,
    required=True,
    help="File path for the list of predicted genes from Prodigal."
)
parser.add_argument(
    "--amr",
    type=str,
    required=True,
    help="File path for the output of AMRFinder."
)
parser.add_argument(
    "--rgi",
    type=str,
    required=True,
    help="File path for the output of RGI."
)
parser.add_argument(
    "--bakta",
    type=str,
    required=True,
    help="File path for the output of Bakta."
)
parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="File path for the output file for this script."
)

# %% Variable Names
# Try to parse args, set defaults if running interactively
try:
    args = parser.parse_args()

    # Name variables from command-line arguments
    genefile = args.genes
    amrfile = args.amr
    rgifile = args.rgi
    baktafile = args.bakta
    outputfile = args.output
except SystemExit:
    # Set default values for interactive development
    print("Running in interactive mode - using default values.")
    genefile = "/Users/aliciasorgen/scratch/BMT101D-7_geneslist.tsv"
    amrfile = "/Users/aliciasorgen/scratch/BMT101D-7.amrfinder.tsv"
    rgifile = "/Users/aliciasorgen/scratch/BMT101D-7.rgi.txt"
    baktafile = "/Users/aliciasorgen/scratch/BMT101D-7.tsv"
    outputfile = "/Users/aliciasorgen/scratch/BMT101D_gene_annotations.tsv"

# %% Defining function to skip rows starting with '#' in Bakta file --------------------------------
def get_skiprows(file_path, to_skip):
    skip_rows = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith(to_skip):
                skip_rows.append(i)
    return skip_rows



# %% Read in gene list data ------------------------------------------------------------------------
genes = pd.read_csv(genefile, header=None, names=["gene_id"], sep='\t')


# %% Process RGI data ------------------------------------------------------------------------------
rgi = pd.read_csv(rgifile, sep='\t')

# Create gene_id column keeping only the first substring from ORF_ID
rgi['gene_id'] = rgi['ORF_ID'].str.split(' ', expand=True)[0]

# Filter rgi to only keep 'Best_Identities' >= 70% & 
# 'Percentage Length of Reference Sequence' >= 50%
rgi = rgi[(rgi['Best_Identities'] >= 70) & 
          (rgi['Percentage Length of Reference Sequence'] >= 50)]

# Keep only the selected columns
rgi = rgi[['gene_id','Model_type','Best_Hit_ARO', 'Drug Class','AMR Gene Family','ARO']]

# Rename the selected columns
rgi = rgi.rename(columns={'Best_Hit_ARO': 'gene_name', 
                          'Drug Class': 'gene_class', 
                          'AMR Gene Family':'gene_product',
                          'ARO':'gene_accession'})

# %% Process AMRFinder data ------------------------------------------------------------------------
amrfinder = pd.read_csv(amrfile, sep='\t')

# Filter amrfinder to only keep '% Identity to reference sequence' >= 70% & 
#   '% Coverage to reference sequence' >= 50%

amrfinder = amrfinder[(amrfinder['% Identity to reference'] >= 70) &
                      (amrfinder['% Coverage of reference'] >= 50)]

# Keep only the selected columns
amrfinder = amrfinder[['Protein id', 'Element symbol', 'Subclass',
                       'Element name','Closest reference accession']]

amrfinder = amrfinder.rename(columns={'Protein id': 'gene_id',
                                      'Element symbol': 'gene_name',
                                      'Subclass': 'gene_class',
                                      'Element name': 'gene_product',
                                      'Closest reference accession': 'gene_accession'})

# %% Process Bakta data ----------------------------------------------------------------------------
skip_rows = get_skiprows(baktafile, "#")
bakta = pd.read_csv(baktafile, sep='\t', skiprows=skip_rows)

bakta = bakta[['ID', 'Gene', 'Product','RefSeq']]
bakta =bakta[bakta['Product'] != "hypothetical protein"]
bakta = bakta.rename(columns={'ID': 'gene_id', 
                              'Gene': 'gene_name', 
                              'Product': 'gene_product',
                              'RefSeq': 'gene_accession'})

# %% Match annotations -----------------------------------------------------------------------------
amr_match = amrfinder.set_index('gene_id').reindex(genes['gene_id'])['gene_name'].values
rgi_match = rgi.set_index('gene_id').reindex(genes['gene_id'])['gene_name'].values
bakta_match = bakta.set_index('gene_id').reindex(genes['gene_id'])['gene_name'].values

# Match AMR class
amr_class = amrfinder.set_index('gene_id').reindex(genes['gene_id'])['gene_class'].values
rgi_class = rgi.set_index('gene_id').reindex(genes['gene_id'])['gene_class'].values

# Match gene family
amr_product = amrfinder.set_index('gene_id').reindex(genes['gene_id'])['gene_product'].values
rgi_product = rgi.set_index('gene_id').reindex(genes['gene_id'])['gene_product'].values
bakta_product = bakta.set_index('gene_id').reindex(genes['gene_id'])['gene_product'].values

# Match gene accession
amr_accession = amrfinder.set_index('gene_id').reindex(genes['gene_id'])['gene_accession'].values
rgi_accession = rgi.set_index('gene_id').reindex(genes['gene_id'])['gene_accession'].values
bakta_accession = bakta.set_index('gene_id').reindex(genes['gene_id'])['gene_accession'].values

#Obtain RGI models
rgi_model = rgi.set_index('gene_id').reindex(genes['gene_id'])['Model_type'].values

# Create the final dataframe
genes['gene_annotation'] = np.where(~pd.isna(rgi_match), rgi_match,
	np.where(~pd.isna(amr_match),amr_match,
	np.where(~pd.isna(bakta_match),bakta_match,"")))
genes['resistance_class'] = np.where(~pd.isna(rgi_match), rgi_class,
	np.where(~pd.isna(amr_match),amr_class,""))
genes['gene_family'] = np.where(~pd.isna(rgi_match), rgi_product, 
	np.where(~pd.isna(amr_match),amr_product,
	np.where(~pd.isna(bakta_match),bakta_product,"")))
genes['gene_accession'] = np.where(~pd.isna(rgi_match), rgi_accession,
	np.where(~pd.isna(amr_match), amr_accession,
	np.where(~pd.isna(bakta_match),bakta_accession,"")))
genes['annotated_with'] = np.where(~pd.isna(rgi_match),rgi_model, 
	np.where(~pd.isna(amr_match),"AMRFinder",
	np.where(~pd.isna(bakta_match),"Bakta","")))

# Write the output to a file
genes.to_csv(outputfile, sep='\t', index=False)