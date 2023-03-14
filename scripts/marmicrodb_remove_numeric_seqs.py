#!/usr/bin/env python3

"""


	# load in FASTA file.

	# Go through records

	# If a sequence starts with numerics- throw it out!

	# only print out good sequences

"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="input MarMicroDB FASTA", type=str)
parser.add_argument("-o", "--output", help="Output FASTA", type=str)
args = parser.parse_args()

# open files
output_fasta = open(args.output, 'w')

# open input_fasta:
input_fasta = open(args.fasta, 'r')

# iterate through the input sequence file:
# get taxid from seq_record.id RefDict and ref_id
for seq_record in SeqIO.parse(input_fasta, "fasta"):
	# if the sequence doesn't begin with a numeric, then pass it on:
	if str(seq_record.seq)[0].isdigit() == False:
		#print("This does not start with a digit:")
		SeqIO.write(seq_record, output_fasta, "fasta")
	else:
		print("This starts with a digit:")
		print(str(seq_record.id))
		print(str(seq_record.seq))
