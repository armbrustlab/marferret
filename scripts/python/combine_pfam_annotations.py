#!/usr/bin/env python3
# AUTHOR: Ryan Groussman

import pandas as pd
import argparse
import os


def handle_arguments():
    description = '''This script aggregates the best-matching Pfam annotations
	in annotations.csv files created by best_pfam.py. Inputs a list of taxIDs to iterate
	through and add the taxID column as a new file into a single output file

        Example usage: ./combine_pfam_annotations.py -f taxIDs.txt MarFERReT.v1.Pfam34.annotations.csv
        '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input', type=str, help='Input list of taxIDs')
    parser.add_argument('output', type=str, help='Output CSV file')
    return parser

def main():
	# parse arguments
	parser = handle_arguments()
	args = parser.parse_args()

	# Read the list of unique identifiers (taxIDs) from the input file
	with open(args.input, 'r') as f:
	    taxIDs = [line.strip() for line in f.readlines()]

	# Initialize an empty DataFrame to store the final output
	output_df = pd.DataFrame(columns=['taxid', 'aa_id', 'pfam_name', 'pfam_id'])

	# Iterate through the list of taxIDs
	for taxID in taxIDs:
	    # Read the CSV file corresponding to the current taxID
	    input_csv = f"{taxID}.Pfam34.annotations.csv"
	    input_df = pd.read_csv(input_csv)

	    # Keep only the required columns
	    input_df = input_df[['aa_id', 'pfam_name', 'pfam_id']]

	    # Add a new column called 'taxid' and populate with the taxID value
	    input_df['taxid'] = taxID

	    # Append the processed input_df to the output_df
	    output_df = pd.concat([output_df, input_df], ignore_index=True)

	# Output the final DataFrame to a CSV file
	output_df.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
