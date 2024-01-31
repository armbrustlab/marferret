#!/usr/bin/env python3

# AUTHOR: Ryan Groussman, PhD

"""
# 1. Load the names of kept ids from KEEP_ID and make set
# 2. Iterate through the MARMICRODB_FASTA
# 3. Parse the mmdb_defline to get the get mmdb_genome
		# Continue only if in keep list
# 4. Inserts a uid for the fasta file as the new defline
# 5. write out a new uid2tax file
"""

def iterate_through_input_fasta(kept_ids, input_fasta):
	"""
	We need to go from this: >150NLHA-308_2182671
	To this: >uidNNNNNNNN 150NLHA-308_2182671
	"""

	# set the sequence counter to zero:
	global count_tracker
	count_tracker = 0

	# iterate through the input sequence file:
	for seq_record in SeqIO.parse(input_fasta, "fasta"):
		"""
		# Defline format for MARMICRODB is >($1)-($2)_$3
		# ex: >150NLHA-308_2182671
		# $1 = 'genome' id
		# $2 = seq_n
		# $3 = NCBI tax_id

		# Check if genome id in kept_ids

		# Output will be:
		# >uid$n ($1)-($2)_$3
		"""
		# get tax_id and genome_id:
		genome_id, tax_id = parse_header(seq_record.id)
		if genome_id in kept_ids:
			# increment the count_tracker:
			count_tracker += 1
			uid = "mmdb" + str(count_tracker)

			# write out sequence with new defline (">mmdb[seq_i] old_defline")
			out_record = SeqRecord(seq_record.seq, id= uid, description=seq_record.id)
			SeqIO.write(out_record, output_fasta, "fasta")

			# write out uid and tax_id to uid2taxid file:
			uid2taxid_line = "NA" + "\t" + uid + "\t" + tax_id + "\t" + "NA\n"
			uid2taxid.write(uid2taxid_line)

def parse_header(seq_record_id):
	"""
	# Defline format for MARMICRODB is >($1)-($2)_$3
	# ex: >150NLHA-308_2182671
	# $1 = 'genome' id
	# $2 = seq_n
	# $3 = NCBI tax_id
	"""
	record_elts = seq_record_id.strip().split("_")
	genome_elts =  record_elts[0].strip(">").split("-")
	seq_n = genome_elts.pop()
	genome_id = "-".join(genome_elts)
	tax_id = str(record_elts[1])
	return genome_id, tax_id


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="input FASTA", type=str)
parser.add_argument("-o", "--output", help="Output FASTA", type=str)
parser.add_argument("-k", "--kept_ids", help="txt file with list of MARMICRODB genome ids to keep", type=str)
parser.add_argument("-t", "--uid2taxid", help="Mapping file of uniq_id to tax_id, tab separated", type=str)

args = parser.parse_args()

# open files and add some headers:
uid2taxid = open(args.uid2taxid, 'w')
output_fasta = open(args.output, 'w')
# We dont need to use this function to add a header as we combine with another file:
# init_uid2taxid_header()

# make a list of the tax_id to keep:
kept_ids = set([str(x).strip() for x in open(args.kept_ids, 'r').readlines()])
print(len(kept_ids))

# open input fasta:
input_fasta = open(args.fasta, 'r')

# count_tracker for counting sequences and assigning uids:
global count_tracker
count_tracker = 0

# perform the main functions on the input fasta
# generating a uid-modified output fasta and uid2tax file
iterate_through_input_fasta(kept_ids, input_fasta)
