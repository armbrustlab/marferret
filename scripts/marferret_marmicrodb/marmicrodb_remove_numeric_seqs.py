#!/usr/bin/env python

# AUTHOR: Ryan Groussman

"""
The purpose of this script is to load in a FASTA file and iterate through
the records, filtering out sequences that start with numeric values and 
printing the passing sequences to an output FASTA file. 

"""

# Import dependences
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Argument parsing for file inputs
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="input MarMicroDB FASTA", type=str)
parser.add_argument("-o", "--output", help="Output FASTA", type=str)
args = parser.parse_args()

# Open input and output FASTA files:
input_fasta = open(args.fasta, 'r')
output_fasta = open(args.output, 'w')

# Iterate through the input sequence file:
for seq_record in SeqIO.parse(input_fasta, "fasta"):
	# if the sequence doesn't begin with a numeric, then pass it on:
	if str(seq_record.seq)[0].isdigit() == False:
		# Write the passing sequence to the out file
		SeqIO.write(seq_record, output_fasta, "fasta")
	else:
		# print the ID and sequence of records with numerics in the sequence field
		print("This starts with a digit:")
		print(str(seq_record.id))
		print(str(seq_record.seq))