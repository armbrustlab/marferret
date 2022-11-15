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
    # initialize counter and mappings
    aa_counter = 0
    uid_mappings = []
    # loop through each taxid
    for i, taxid in enumerate(meta_df['tax_id'].unique()):
        # subset metadata
        subset_df = meta_df[meta_df['tax_id'] == taxid]
        print('Combining {} reference(s) of taxid {}'.format(len(subset_df), taxid))
        # open new outfile
        with open('{}/{}.combined.faa'.format(args.output_dir, taxid), 'w') as outfile:
            # iterate through input fasta files
            for input_fasta in subset_df['aa_fasta']:
                # iterate through SeqRecords
                for seq_record in SeqIO.parse('{}/{}'.format(args.input_dir, input_fasta), 'fasta'):
                    # make new uid for sequence
                    aa_counter += 1 
                    uid = 'mft{:0>10}'.format(aa_counter)
                    uid = 'mft{:d}'.format(aa_counter)
                    # write renamed SeqRecord to outfile
                    out_record = SeqRecord(seq_record.seq, id=uid, description='')
                    SeqIO.write(out_record, outfile, 'fasta')
                    # record uid mappings
                    uid_mappings.append(
                        {
                            'marferret_id': uid, 
                            'source_defline': seq_record.id.strip(), 
                            'taxid': taxid
                        }
                    )
    # save uid mapping files
    print('Assembling reference mapping files')
    mapping_df = pd.DataFrame(uid_mappings)
    # save uid2def.csv file
    print('Saving reference mapping files')
    mapping_df[['marferret_id', 'source_defline']].to_csv(
        '{}/uid2def.csv'.format(args.output_dir), 
        index=False
    )
    # save uid2tax.tab file
    mapping_df['accession'] = 'NA'
    mapping_df['gi'] = 'NA'
    mapping_df.rename(columns={'marferret_id': 'accession.version'}, inplace=True)
    mapping_df[['accession', 'accession.version', 'taxid', 'gi']].to_csv(
        '{}/uid2tax.tab'.format(args.output_dir), 
        sep='\t', 
        index=False
    )
    print('Finished!')
    
if __name__ == "__main__":
    main()
    