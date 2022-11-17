#!/usr/bin/env python3

# imports
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def handle_arguments():
    description = '''This script groups the input fasta files by tax_id, 
    outputting a single fasta file for each tax_id group. It also renames 
    each of the SeqRecords with a unique identifier, and outputs mapping 
    files that match each new unique id to its old id and tax_id. 
    
    Example usage: ./group_by_taxid.py fasta/dir/in metadata.csv -o fasta/dir/out
    '''
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        'input_dir', 
        type=str, 
        help='Directory containing input fasta sequences.'
    )
    parser.add_argument(
        'metadata', 
        type=str, 
        help='Metadata csv file containing "aa_fasta" and "tax_id" fields.'
    )
    parser.add_argument(
        '-o', '--output_dir', 
        type=str, 
        help='Directory to which output fasta and reference files will be written.'
    )
    return parser

def main():
    # parse arguments
    parser = handle_arguments()
    args = parser.parse_args()
    # read in metadata
    meta_df = pd.read_csv(args.metadata)
    # open and initialize mapping files
    uid2tax = open('{}/uid2tax.tab'.format(args.output_dir), 'w')
    uid2tax.write('accession\taccession.version\ttaxid\tgi\n')
    uid2def = open('{}/uid2def.csv'.format(args.output_dir), 'w')
    uid2def.write('uniq_id,source_defline\n')
    # initialize counter
    aa_counter = 0
    # loop through each taxid
    for i, taxid in enumerate(meta_df['tax_id'].unique()):
        # subset metadata
        subset_df = meta_df[meta_df['tax_id'] == taxid]
        print('Combining {} reference(s) of taxid {}'.format(len(subset_df), taxid), flush=True)
        # open new outfile
        with open('{}/{}.combined.faa'.format(args.output_dir, taxid), 'w') as outfile:
            # iterate through input fasta files
            for input_fasta in subset_df['aa_fasta']:
                # iterate through SeqRecords
                for seq_record in SeqIO.parse('{}/{}'.format(args.input_dir, input_fasta), 'fasta'):
                    # make new uid for sequence
                    aa_counter += 1 
                    uid = 'mft{:d}'.format(aa_counter)
                    # write renamed SeqRecord to outfile
                    out_record = SeqRecord(seq_record.seq, id=uid, description='')
                    SeqIO.write(out_record, outfile, 'fasta')
                    # write uid mappings to mapping files
                    uid2tax.write('NA\t{}\t{}\tNA\n'.format(uid, taxid))
                    uid2def.write('{},{}\n'.format(uid, seq_record.id.strip()))
    # close file handlers
    uid2tax.close()
    uid2def.close()
    print('Finished!', flush=True)
    
if __name__ == "__main__":
    main()
    