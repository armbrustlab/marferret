#!/usr/bin/env python3


# uniq_id_and_mapping.py

"""
Input:

	REF_DAT="/mnt/nfs/projects/ryan/EukRefDB_2021/metadata/EukRefDB_current.aa_paths.csv"

	Colnames and first entry:
	ref_id	Genus_species	6tr_dir	6tr	NCBI tax_id
	1	Minidiscus_spinulatus	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-01_Minidiscusspinulatus-Trinity.6tr.fasta	2593073

	FASTA files with path & file name linked in REF_DAT
	FASTA files will be loaded by the 'taxID' variable.
	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-06_RCC4590-Minidiscussp-Trinity.6tr.fasta
	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-07_RCC4582-Minidiscussp-Trinity.6tr.fasta
	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-10_RCC4584-Minidiscussp-Trinity.6tr.fasta

Process:
1. Load in REF_DAT init_ref_dat_dictionary(ref_dat)
	1. Create a dictionary with the taxID as keys and
		containing dictionary with:
		 ref_id:
		 	ref_id aa_fasta tax_id entries

2. Use the dictionary to iterate through each taxID.
	For each taxID,
		For each ref_id in each taxID:
			Load the linked fasta file,
			Create output file ($TAXID.fasta) in combined folder
			Reset sequence counter to 0
				For each sequence in each ref_id:
					Advance the sequence counter

1. Parse the defline for each sequence
	a. It will get a unique identifier (uid#) based on an incremental count_tracker
	b. taxID will be gathered from [taxID][ref_id][tax_id]
	c. Create uid: mft[iterator]

2. Outputs:
	a. New fasta file (taxID.fasta )with uid# in defline: >mft#### defline
		# scrap the defline to save space actually; it won't be used unless
		# necessary. some of them are very large. users can refer to Mapping
		# file 2 to get the link between uniq_id and original defline
	b. Mapping file 1: uniq_id to tax_id (tab separated, for DIAMOND input)
		()
	c. Mapping file 2: uniq_id to full defline
		# the MarFERReT uniq_id and the original source defline
"""

# Output directory:

def init_ref_dat_dictionary(ref_dat):
	# 1. Load in REF_DAT

	# 	1. Create a dictionary with the taxID as keys and
	# 		containing dictionary with:
	# 		 ref_id:
	# 		 ref_id 6tr_dir	6tr	NCBI tax_id entries

	# initiate empty dictionary:
	RefDict = {}
	for row in ref_dat:
		ref_id = row['ref_id']
		aa_fasta = row['aa_fasta']
		tax_id = row['tax_id']

		# if taxID is not a key in dict, create it as a key linking to a dict:
		if tax_id not in RefDict.keys():
			RefDict[tax_id] = {}
		# create ref_id as entry under the tax_id, create it as a dict:
		RefDict[tax_id][ref_id] = {}
		# populate RefDict[tax_id][ref_id] with values
		RefDict[tax_id][ref_id]['tax_id'] = tax_id
		RefDict[tax_id][ref_id]['path'] = "/".join([fasta_dir,aa_fasta])

	return RefDict


def create_taxid_file(RefDict,tax_id):

	# take in RefDict and a tax_id entry.

	# Create an output fasta specifically for the tax_id:
	taxid_fasta_path = ".".join([tax_id,"combined.faa"])

	global output_fasta
	output_fasta = open(taxid_fasta_path, 'w')

	# iterate_through_fasta(tax_id, input_fasta)
	for ref_id in RefDict[tax_id].keys():
		print("Opening reference # " + ref_id)
		print(RefDict[tax_id][ref_id]['path'])
		iterate_through_input_fasta(RefDict, tax_id, ref_id)

def iterate_through_input_fasta(RefDict, tax_id, ref_id):

	# open input_fasta:
	input_fasta_path = RefDict[tax_id][ref_id]['path']
	input_fasta = open(input_fasta_path, 'r')

	## set the sequence counter to zero:
	global count_tracker
	#count_tracker = 0

	# iterate through the input sequence file:
	for seq_record in SeqIO.parse(input_fasta, "fasta"):

		# increment the count_tracker:
		count_tracker += 1
		uid = "mft" + str(count_tracker)

		# write out sequence with new defline (">mft[seq_i]")
		out_record = SeqRecord(seq_record.seq, id=uid, description="")
		SeqIO.write(out_record, output_fasta, "fasta")
		# write out uid and tax_id to uid2taxid file:
		uid2taxid_line = "NA" + "\t" + uid + "\t" + tax_id + "\t" + "NA\n"
		uid2taxid.write(uid2taxid_line)
		# write out uid and old_defline to uid2defline file:
		old_defline = seq_record.id.strip()
		uid2defline_line = uid + "," + old_defline + "\n"
		uid2defline.write(uid2defline_line)


def init_uid2taxid():

	# write out the header:
	outheader = "accession\taccession.version\ttaxid\tgi\n"
	uid2taxid.write(outheader)

def init_uid2defline():

	# write out the header:
	outheader = "uniq_id,source_defline\n"
	uid2defline.write(outheader)

import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--uid2taxid", help="Mapping file of uniq_id to tax_id, tab separated", type=str)
parser.add_argument("-r", "--ref_dat", help="CSV containing ref_id,tax_id,aa_path", type=str)
parser.add_argument("-c", "--uid2defline", help="Mapping file of uniq_id to original defline, comma separated", type=str)
parser.add_argument("-d", "--fasta_dir", help="Path of input FASTA files", type=str)
args = parser.parse_args()

# count_tracker for counting sequences and assigning uids:
global count_tracker
count_tracker = 0

# load in the refdb data csv
print("Loading and parsing reference db data...")
fasta_dir = args.fasta_dir
ref_dat = csv.DictReader(open(args.ref_dat))
RefDict = init_ref_dat_dictionary(ref_dat)

# open files and add some headers:
uid2taxid = open(args.uid2taxid, 'w')
init_uid2taxid()
uid2defline = open(args.uid2defline, 'w')
init_uid2defline()

# now iterate through each tax_id in RefDict and create a fasta for each:
for tax_id in RefDict.keys():
	create_taxid_file(RefDict, tax_id)

print("Finished!")
