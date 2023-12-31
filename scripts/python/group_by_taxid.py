#!/usr/bin/env python3

# imports
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def handle_arguments():
    description = '''This script groups the input fasta files by tax_id, 
    outputting a single fasta file for each tax_id group. It also renames 
    each of the SeqRecords with a unique amino acid identifier `aa_id`, 
    and outputs mapping files that match each new aa_id to its corresponding 
    entry_id, tax_id, as well as the name of the sequence in the original input
    fasta file. 
    
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
        help='Metadata csv file containing "entry_id", "marferret_name" and "tax_id" fields.'
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
    aa_tax = open('{}/taxonomies.tab'.format(args.output_dir), 'w')
    aa_tax.write('accession\taccession.version\ttaxid\tgi\n')
    aa_info = open('{}/proteins_info.tab'.format(args.output_dir), 'w')
    aa_info.write('aa_id\tentry_id\tsource_defline\n')
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
            for entry in subset_df.iterrows():
                # get important entry information
                marferret_name = entry[1]['marferret_name']
                entry_id = entry[1]['entry_id']
                input_fasta = str(entry_id) + "_" + marferret_name + ".faa"
                # iterate through SeqRecords
                for seq_record in SeqIO.parse('{}/{}'.format(args.input_dir, input_fasta), 'fasta'):
                    # make new uid for sequence
                    aa_counter += 1 
                    uid = 'mft{:0>10}'.format(aa_counter)
                    # write renamed SeqRecord to outfile
                    out_record = SeqRecord(seq_record.seq, id=uid, description='')
                    SeqIO.write(out_record, outfile, 'fasta')
                    # write uid mappings to mapping files
                    aa_tax.write('NA\t{}\t{}\tNA\n'.format(uid, taxid))
                    aa_info.write('{}\t{}\t{}\n'.format(uid, entry_id, seq_record.id.strip()))
    # close file handlers
    aa_tax.close()
    aa_info.close()
    print('Finished!', flush=True)
    
if __name__ == "__main__":
    main()
    
