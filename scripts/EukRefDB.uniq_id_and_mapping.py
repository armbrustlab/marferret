#!/usr/bin/env python3


# EukRefDB.uniq_id_and_mapping.py

"""
Input:

	REF_DAT="/mnt/nfs/projects/ryan/EukRefDB_2021/metadata/EukRefDB_current.aa_paths.csv"

	Colnames and first entry:
	ref_id	Genus_species	6tr_dir	6tr	NCBI tax_id
	1	Minidiscus_spinulatus	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-01_Minidiscusspinulatus-Trinity.6tr.fasta	2593073

	FASTA files with path & file name linked in REF_DAT
	FASTA files will be loaded by the 'Genus_species' variable.
	e.g.: for Genus_species = Minidiscus_sp
	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-06_RCC4590-Minidiscussp-Trinity.6tr.fasta
	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-07_RCC4582-Minidiscussp-Trinity.6tr.fasta
	/mnt/nfs/projects/ryan/EukRefDB_2021/source/Guijardo_2021/6tr	ThTSP-10_RCC4584-Minidiscussp-Trinity.6tr.fasta

Process:
1. Load in REF_DAT init_ref_dat_dictionary(ref_dat)
	1. Create a dictionary with the Genus_species as keys and
		containing dictionary with:
		 ref_id:
		 ref_id 6tr_dir	6tr	NCBI tax_id entries

2. Use the dictionary to iterate through each Genus_species.
	For each Genus_species,
		For each ref_id in each Genus_species:
			Load the linked fasta file,
			Create output file (Genus_species.fasta) in combined folder
			Reset sequence counter to 0
				For each sequence in each ref_id:
					Advance the sequence counter

1. Parse the defline for each sequence
	a. It will get a unique identifier (uid#) based on an incremental count_tracker and ref_id
	b. taxid will be gathered from [Genus_species][ref_id][tax_id]
	c. Create uid: ref[ref_id]_seq[iterator]

2. Outputs:
	a. New fasta file (Genus_species.fasta )with uid# in defline: >uid#### defline
	b. Mapping file 1: uniq_id to tax_id (tab separated, for DIAMOND input)
	c. Mapping file 2: uniq_id to full defline

"""

# Output directory:

def init_ref_dat_dictionary(ref_dat):
	# 1. Load in REF_DAT

	# 	1. Create a dictionary with the Genus_species as keys and
	# 		containing dictionary with:
	# 		 ref_id:
	# 		 ref_id 6tr_dir	6tr	NCBI tax_id entries

	# initiate empty dictionary:
	RefDict = {}
	for row in ref_dat:
		ref_id = row['ref_id']
		Genus_species = row['Genus_species']
		fasta_dir = row['6tr_dir']
		fasta_file = row['6tr']
		tax_id = row['tax_id']

		# if Genus_species is not a key in dict, create it as a key linking to a dict:
		if Genus_species not in RefDict.keys():
			RefDict[Genus_species] = {}
		# if the tax_id is not an entry under the Genus_species, create it as a dict:
		if ref_id not in RefDict.keys():
			RefDict[Genus_species][ref_id] = {}
		# populate RefDict[Genus_species][ref_id] with values
		RefDict[Genus_species][ref_id]['tax_id'] = tax_id
		RefDict[Genus_species][ref_id]['path'] = "/".join([fasta_dir,fasta_file])

	return RefDict


def create_species_file(RefDict,Genus_species):

	# take in RefDict and a Genus_species entry.

	# Create an output fasta specifically for the Genus_species:
	species_fasta_path = ".".join([Genus_species,"combined.aa.fasta"])

	global output_fasta
	output_fasta = open(species_fasta_path, 'w')

	# iterate_through_fasta(Genus_species, input_fasta)
	for ref_id in RefDict[Genus_species].keys():
		print("Opening reference # " + ref_id)
		print(RefDict[Genus_species][ref_id]['path'])
		iterate_through_input_fasta(RefDict, Genus_species, ref_id)

def iterate_through_input_fasta(RefDict, Genus_species, ref_id):

	# open input_fasta:
	input_fasta_path = RefDict[Genus_species][ref_id]['path']
	input_fasta = open(input_fasta_path, 'r')

	# set the sequence counter to zero:
	global count_tracker

	count_tracker = 0

	# iterate through the input sequence file:
	for seq_record in SeqIO.parse(input_fasta, "fasta"):
		# get uid, taxid, and old_defline from seq_record.id RefDict and ref_id
		uid, tax_id, old_defline = parse_header(seq_record.id, RefDict, Genus_species, ref_id)

		# write out sequence with new defline (">ref[ref_id]_seq[seq_i] old_defline")
		out_record = SeqRecord(seq_record.seq, id= uid, description=old_defline)
		SeqIO.write(out_record, output_fasta, "fasta")
		# write out uid and tax_id to uid2taxid file:
		uid2taxid_line = "NA" + "\t" + uid + "\t" + tax_id + "\t" + "NA\n"
		uid2taxid.write(uid2taxid_line)
		# write out uid and old_defline to uid2defline file:
		uid2defline_line = uid + "," + old_defline + "\n"
		uid2defline.write(uid2defline_line)

def parse_header(seq_record_id, RefDict, Genus_species, ref_id):

	# get uid, taxid, and defline from seq_record.id RefDict and ref_id
	# we will create the defline uid in this format:
	# ref[ref_id]_seq[seq_i]

	# increment the count_tracker:
	global count_tracker
	count_tracker += 1
	uid = "ref" + str(ref_id) + "_seq" + str(count_tracker)
	# split by underscore and grab the taxid:
	old_defline = seq_record_id.strip()

	tax_id = RefDict[Genus_species][ref_id]['tax_id']

	return uid, tax_id, old_defline

def init_uid2taxid():

	# write out the header:
	outheader = "accession\taccession.version\ttaxid\tgi\n"
	uid2taxid.write(outheader)

def init_uid2defline():

	# write out the header:
	outheader = "uniq_id,full_defline\n"
	uid2defline.write(outheader)

import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--uid2taxid", help="Mapping file of uniq_id to tax_id, tab separated", type=str)
parser.add_argument("-r", "--ref_dat", help="CSV containing ref_id,Genus_species,6tr_dir,6tr,NCBI tax_id", type=str)
parser.add_argument("-c", "--uid2defline", help="Mapping file of uniq_id to original defline, comma separated", type=str)
args = parser.parse_args()

# count_tracker for counting sequences and assigning uids:
global count_tracker
count_tracker = 0

# load in the refdb data csv
print("Loading and parsing reference db data...")
ref_dat = csv.DictReader(open(args.ref_dat))
RefDict = init_ref_dat_dictionary(ref_dat)

# open files and add some headers:
uid2taxid = open(args.uid2taxid, 'w')
init_uid2taxid()
uid2defline = open(args.uid2defline, 'w')
init_uid2defline()

# now iterate through each Genus_species in RefDict and create a fasta for each:
for Genus_species in RefDict.keys():
	create_species_file(RefDict, Genus_species)

print("Finished!")
