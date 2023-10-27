#!/usr/bin/env python3

# Ryan Groussman, PhD
# October 2023
# Armbrust lab
# School of Oceanography
# University of Washington

import os
import pandas as pd


def collect_input_data(this_entry):
	# gather and process the input files
		
	# input the hmm file (seq to pfam) 
	hmmfile_path = "../../headers/" + this_entry + "_ribo_fasta_headers.csv"
	hmmfile = pd.read_csv(hmmfile_path, names=hmm_colnames)
	# trim the pfam_id:
	hmmfile['pfam_id'] = hmmfile['pfam_id'].str.split('.').str[0]
	
	# and lca file (seq to taxid)
	lcafile_path = this_entry + "_ribosomal.pfam_full.lca.tab"
	lcafile = pd.read_csv(lcafile_path, sep='\t', names=lca_colnames)
	
	# merge to one seq*taxid*pfam file
	merge_dat = pd.merge(hmmfile[["aa_id","pfam_name","pfam_id"]], lcafile[["aa_id", "tax_id"]], on="aa_id",how="left")
	
	# filter to the 41 core RPs
	merge_dat_ctg = merge_dat[merge_dat['pfam_id'].isin(ctg_pfam_ids)]
	
	# return merge_dat
	return merge_dat

def init_entry_dict(this_entry):
	# This function builds a dictionary for predefined taxonomic lineages
	# using the NCBI taxID as key to a dictionary with the taxa name, rank, and a
	# parameter 'count' for number of sequences
	entry_id = int(this_entry.split('_')[0])

	EntryDict = {'entry_handle' : this_entry, 'entry_id': entry_id,
	'tax_bin' : {
	554915 : {'tax_name' : 'Amoebozoa', 'rank' : 'superkingdom_', 'count':0},
	28009: {'tax_name' : 'Choanoflagellata', 'rank' : 'class', 'count':0},
	5878: {'tax_name' : 'Ciliophora', 'rank' : 'phylum', 'count':0},
	877183: {'tax_name' : 'Colpodellida', 'rank' : 'superkingdom___', 'count':0},
	3027: {'tax_name' : 'Cryptophyceae', 'rank' : 'class', 'count':0},
	2864: {'tax_name' : 'Dinophyceae', 'rank' : 'class', 'count':0},
	33682: {'tax_name' : 'Euglenozoa', 'rank' : 'phylum', 'count':0},
	4751: {'tax_name' : 'Fungi', 'rank' : 'kingdom', 'count':0},
	38254: {'tax_name' : 'Glaucocystophyceae', 'rank' : 'class', 'count':0},
	2830: {'tax_name' : 'Haptophyta', 'rank' : 'phylum', 'count':0},
	5752: {'tax_name' : 'Heterolobosea', 'rank' : 'phylum', 'count':0},
	759891: {'tax_name' : 'Palpitomonas', 'rank' : 'genus', 'count':0},
	2497438: {'tax_name' : 'Perkinsozoa', 'rank' : 'phylum', 'count':0},
	543769: {'tax_name' : 'Rhizaria', 'rank' : 'superkingdom__', 'count':0},
	2763: {'tax_name' : 'Rhodophyta', 'rank' : 'phylum', 'count':0},
	33634: {'tax_name' : 'Stramenopiles', 'rank' : 'superkingdom__', 'count':0},
	33090: {'tax_name' : 'Viridiplantae', 'rank' : 'kingdom', 'count':0},
	2: {'tax_name' : 'Bacteria', 'rank' : 'superkingdom', 'count':0},
	2157: {'tax_name' : 'Archaea', 'rank' : 'superkingdom', 'count':0},
	10239: {'tax_name' : 'Viruses', 'rank' : 'superkingdom', 'count':0},
	0 : {'tax_name' : 'Other', 'rank' : 'NA', 'count':0},
	}, 
	'ranks' : {
	'superkingdom' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'kingdom' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'phylum' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'class' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'order' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'family' :{'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'genus' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'species' :{'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0}
	}}
	return EntryDict

def get_expected_ranks(this_entry):
	# get the entry_id from this_entry:
	entry_id = int(this_entry.split('_')[0])
	# use metadata.csv to get the entry_taxid
	entry_taxid = metadat.loc[metadat['candidate_id']==entry_id, 'tax_id'].item()
	# add tax_id to EntryDict
	EntryDict['tax_id'] = int(entry_taxid)
	
	# get taxa.csv dat from the marferret taxa.csv file:
	entry_dat = mft_taxa_csv[mft_taxa_csv['tax_id'] == entry_taxid]
	
	# set the ranks in EntryDict
	for rank in EntryDict['ranks'].keys():
		rank_taxid = entry_dat.loc[entry_dat['tax_id']==entry_taxid, rank].item()
		try: rank_taxid = int(rank_taxid)
		except ValueError: rank_taxid = "NA"
		if rank_taxid != "NA":
			EntryDict['ranks'][rank]['expected'] = rank_taxid
	
	# set expected group from the test_ranks:
	EntryDict['expected_id'] = "None"
	EntryDict['expected_name'] = "None"
	for test_rank in test_ranks:
		test_id = entry_dat.loc[entry_dat['tax_id']==entry_taxid, test_rank].item()
		try: test_id = int(test_id)
		except ValueError: test_id = "NA"
		for tax_bin in EntryDict['tax_bin'].keys():
			if tax_bin == test_id:
				EntryDict['expected_id'] = test_id
				EntryDict['expected_name'] = EntryDict['tax_bin'][tax_bin]['tax_name']

def check_if_defined_taxa(query_taxid):
	# if tax_id is '0' then assign Unknown values and return
	if query_taxid == 0:
		EntryDict['tax_bin'][0]['count'] += 1
		return ""
	# get the data for this tax_id from the taxa.csv file
	taxon_dat = taxa_csv[taxa_csv['tax_id'] == query_taxid]
	# iterate thru ranks:
	for test_rank in test_ranks:
		test_id = taxon_dat.loc[taxon_dat['tax_id']==query_taxid, test_rank].item()
		try: test_id = int(test_id)
		except ValueError: test_id = "NA"
		if test_id in EntryDict['tax_bin'].keys():
			# increment counter:
			EntryDict['tax_bin'][test_id]['count'] += 1

def check_expected_ranks(query_taxid):
	# check if query_taxid is in EntryDict expected ranks
	# get query_taxid rank
	query_dat = taxa_csv[taxa_csv['tax_id'] == query_taxid]
	# check if query rank ids match expected rank ids:
	for rank in EntryDict['ranks'].keys():
		try: rank_taxid = query_dat.loc[query_dat['tax_id']==query_taxid, rank].item()
		except ValueError: continue
		try: rank_taxid = int(rank_taxid)
		except ValueError: continue
		if rank_taxid == "NA":
			continue
		if rank_taxid == EntryDict['ranks'][rank]['expected']:
			EntryDict['ranks'][rank]['in_rank'] += 1
		elif rank_taxid != EntryDict['ranks'][rank]['expected']:
			EntryDict['ranks'][rank]['not_in_rank'] += 1

def process_entry_pfam(pfam_id, entry_data): 
	# subset entry/pfam df:
	sub_seq_dat =  entry_data[entry_data['pfam_id'] == pfam_id]
	# for each seq in this_pfam
	for query_taxid in sub_seq_dat['tax_id']:
		# get n in each of the set groups (hapto, dino, etc)
		# and update if necessary
		check_if_defined_taxa(query_taxid)
		check_expected_ranks(query_taxid)
	
def entry_stats():
	# n sequences for entry:
	EntryDict['n_seqs'] = len(entry_data['aa_id'])
	EntryDict['n_pfams'] = len(entry_data['pfam_id'].unique())
	expected_bin = EntryDict['expected_id']
	ingroup_count = "NA"
	outgroup_count = "NA"
	contam_pct = "NA"
	if expected_bin != "None":
		ingroup_count = EntryDict['tax_bin'][expected_bin]['count'] 
		outgroup_count = 0
		for tax_bin in EntryDict['tax_bin'].keys():
			if tax_bin != expected_bin:
				outgroup_count += EntryDict['tax_bin'][tax_bin]['count']
		contam_pct = round((outgroup_count/(outgroup_count+ingroup_count)*100),2)
	EntryDict['ingroup_count'] = str(ingroup_count)
	EntryDict['outgroup_count']	= str(outgroup_count)
	EntryDict['contam_pct'] = str(contam_pct)

def init_output_file():
	# open the output file and write header
	outcsv_file = open(outcsv_path,'a')
	# output columns:
	out_cols = ["entry_handle", "entry_id", "tax_id", "n_seqs", "n_pfams", "tax_group", "contam_pct", "Amoebozoa", "Choanoflagellata", "Ciliophora", "Colpodellida", "Cryptophyceae", "Dinophyceae", "Euglenozoa", "Fungi", "Glaucocystophyceae", "Haptophyta", "Heterolobosea", "Palpitomonas", "Perkinsozoa", "Rhizaria", "Rhodophyta", "Stramenopiles", "Viridiplantae", "Bacteria", "Archaea", "Viruses", "Other", 'in_superkingdom', 'in_kingdom', 'in_phylum', 'in_class', 'in_order', 'in_family', 'in_genus', 'in_species']
	
	out_header = ",".join(out_cols) + "\n"
	outcsv_file.write(out_header)
	return outcsv_file

def write_results(this_entry, EntryDict, outcsv_file):
	# write out the results for each entry
	out_line = [this_entry]
	out_line.append(str(EntryDict['entry_id']))
	out_line.append(str(EntryDict['tax_id']))
	out_line.append(str(EntryDict['n_seqs']))
	out_line.append(str(EntryDict['n_pfams']))
	out_line.append(str(EntryDict['expected_name']))
	out_line.append(str(EntryDict['contam_pct']))
	
	for taxbin in EntryDict['tax_bin'].keys():
		out_line.append(str(EntryDict['tax_bin'][taxbin]['count']))
	for rank in EntryDict['ranks'].keys():
		out_line.append(str(EntryDict['ranks'][rank]['in_rank']))
	# write out:
	outline_ready = ",".join(out_line) + "\n"
	outcsv_file.write(outline_ready)


os.chdir("marferret-main/data/pfam/ribosomal_proteins/mft_ribo_sequences/lca_full/")

# output file:
outcsv_path = "marferret-main/data/pfam/ribosomal_proteins/MarFERReT.QC_estimates.csv"

# input taxa.csv
taxa_csv = pd.read_csv("../../Pfam.ribosomal_all_proteins.taxa.csv")

# the smaller mft taxa.csv too
mft_taxa_csv = pd.read_csv("marferret/marferret-main/data/MarFERReT.v1.taxa.csv")

# input metadata.csv
metadat = pd.read_csv("marferret/marferret-main/data/MarFERReT.v1.entry_curation.csv")
metadat.shape # (902, 22)

# hmm.csv column names:
hmm_colnames = ["aa_id", "pfam_name", "pfam_id", "eval", "bitscore"]
lca_colnames = ["aa_id", "tax_id", "eval"]

# list of column/phylum names to use
test_ranks = ["superkingdom_", "superkingdom__", "superkingdom___", "class", "phylum", "kingdom", "genus", "phylum"]

# from [[tech_validation]] ##### Get set of MarFERReT RPs from CTGs
ctg_pfam_ids = ['PF00164', 'PF00177', 'PF00189', 'PF00203', 'PF00237', 'PF00238','PF00252', 'PF00276', 'PF00297', 'PF00312', 'PF00318', 'PF00333', 'PF00338', 'PF00347', 'PF00380', 'PF00400', 'PF00410', 'PF00411','PF00416', 'PF00428', 'PF00466', 'PF00501', 'PF00572', 'PF00573','PF00583', 'PF00687', 'PF00828', 'PF00849', 'PF01015', 'PF01159','PF01189', 'PF01196', 'PF01201', 'PF01246', 'PF01248', 'PF01280','PF01728', 'PF01778', 'PF01926', 'PF03947', 'PF17135']



if __name__ == "__main__":

	# input list of entry handles:
	with open("mft_ribo_sequences.handles.txt") as f:
		entry_handles = f.read().splitlines()

	outcsv_file = init_output_file()

	# for each entry:
	for this_entry in entry_handles:	
		print(this_entry)
		# gather and process the input files
		entry_data = collect_input_data(this_entry)
		# init entry dict:
		EntryDict = init_entry_dict(this_entry)
		# gather expected ranks:
		get_expected_ranks(this_entry)
		# iterate through pfams:
		entry_pfam_list = entry_data['pfam_id'].unique()
		for this_pfam in entry_pfam_list:
			process_entry_pfam(this_pfam, entry_data)
		# run statistics on expected and counted taxids
		entry_stats()
		# then write out EntryDict results:
		write_results(this_entry, EntryDict, outcsv_file)

	# need this to write out output
	outcsv_file.close()
