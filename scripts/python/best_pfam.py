#!/usr/bin/env python
# AUTHOR: Stephen Blaskowski

import pandas as pd
import argparse
import os

def handle_arguments():
    description = '''This script parses the pfam domtblout file, 
        selects the top scoring annotation for each protein listed in the
        file, and then outputs a csv of the parsed and down-selected data.

        Example usage: ./best_pfam.py -f domtblout.tab annotations.csv
        '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('input', type=str, help='Input domtblout file')
    parser.add_argument('output', type=str, help='Output file')
    return parser

def main():
    # parse arguments
    parser = handle_arguments()
    args = parser.parse_args()
    # check output path
    output_path = args.output
    if os.path.isfile(output_path):
        raise FileExistsError('A file by the name of {} already exists.'.format(output_path))
    # parse pfam domtblout data
    print('Parsing domtblout', flush=True)
    df = pd.read_csv(
        args.input, 
        engine='python', 
        sep='\s+', 
        skiprows=[0, 1, 2], 
        skipfooter=10, 
        usecols=[0, 3, 4, 6, 7],
        names=['aa_id', 'pfam_name', 'pfam_id', 'pfam_eval', 'pfam_score']
    )
    # take the best result for each aa_id
    print('Selecting best annotation for each protein', flush=True)
    best_annot_df = df.iloc[df.groupby(['aa_id']).pfam_score.idxmax()]
    best_annot_df.reset_index(drop=True, inplace=True)
    
    # save results to a csv
    best_annot_df.to_csv(args.output, index=False)
    print('Finished')

if __name__ == "__main__":
    main()
