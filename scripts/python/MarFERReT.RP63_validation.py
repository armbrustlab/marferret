#!/usr/bin/env python3

# Ryan Groussman, PhD
# October 2023
# Armbrust lab
# School of Oceanography
# University of Washington

import os
import pandas as pd

# define local marferret directory:
marferret_dir = os.getcwd()

# output file:
outcsv_path = marferret_dir + "MarFERReT.v1.RP63_QC_estimates.csv"

# input taxa.csv for full set of taxIDs from unitprokb
taxa_csv = pd.read_csv(marferret_dir + "/data/pfam/ribosomal_proteins/Pfam.ribosomal_all_proteins.taxa.csv")

# input the taxa.csv for the marferret taxIDs
mft_taxa_csv = pd.read_csv(marferret_dir + "/data/diamond/ncbi/MarFERReT.v1.taxa.csv")

# input entry metadata.csv
metadat = pd.read_csv(marferret_dir + "/data/MarFERReT.v1.metadata.csv")

# hmm.csv column names:
hmm_colnames = ["aa_id", "pfam_name", "pfam_id", "eval", "bitscore"]
lca_colnames = ["aa_id", "tax_id", "eval"]

# list of column/phylum names to use
test_ranks = ["superkingdom_", "superkingdom__", "superkingdom___", "class", "phylum", "kingdom", "genus", "phylum"]

# Set of ribosomal protein Pfam IDs present in >90% of entries:
ctg_pfam_ids63 =["PF01248", "PF01201", "PF00687", "PF00828", "PF00572", "PF00380", "PF00297", "PF00237", "PF00573", "PF00164", "PF03947", "PF00347", "PF01246", "PF00177", "PF01778", "PF00238", "PF00318", "PF01015", "PF17144", "PF00252", "PF00203", "PF16906", "PF00428", "PF00411", "PF01655", "PF01280", "PF00416", "PF00410", "PF17135", "PF00189", "PF01092", "PF01159", "PF00827", "PF01775", "PF01777", "PF01282", "PF00831", "PF00276", "PF01196", "PF00338", "PF01294", "PF01251", "PF00900", "PF01929", "PF01199", "PF00833", "PF00333", "PF01090", "PF01157", "PF00466", "PF03297", "PF01158", "PF00312", "PF01283", "PF01247", "PF01198", "PF08069", "PF00935", "PF00542", "PF03946", "PF00829", "PF01776", "PF00366"]

# filtered entries
# entries removed for LOW_SEQS or LOW_PFAMs
removed_entries = ["195", "262", "426", "730", "255", "267", "276", "280", "282", "302", "321", "371", "383", "394", "425", "446", "472", "560", "629", "665", "669", "706", "712", "818", "847", "870", "901", "902"]


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
	
	# filter to the core RPs
	merge_dat_ctg = merge_dat[merge_dat['pfam_id'].isin(ctg_pfam_ids)]
	
	# return merge_dat
	return merge_dat

def init_entry_dict(this_entry):
	entry_id = int(this_entry.split('_')[0])
	EntryDict = {'entry_handle' : this_entry, 'entry_id': entry_id,
	'tax_bin' : {
	554915 : {'tax_name' : 'Amoebozoa', 'rank' : 'superkingdom_', 'count':0},
	5878: {'tax_name' : 'Ciliophora', 'rank' : 'phylum', 'count':0},
	877183: {'tax_name' : 'Colpodellida', 'rank' : 'superkingdom___', 'count':0},
	3027: {'tax_name' : 'Cryptophyceae', 'rank' : 'class', 'count':0},
	2864: {'tax_name' : 'Dinophyceae', 'rank' : 'class', 'count':0},
	33682: {'tax_name' : 'Euglenozoa', 'rank' : 'phylum', 'count':0},
	#4751: {'tax_name' : 'Fungi', 'rank' : 'kingdom', 'count':0},
	38254: {'tax_name' : 'Glaucocystophyceae', 'rank' : 'class', 'count':0},
	2830: {'tax_name' : 'Haptophyta', 'rank' : 'phylum', 'count':0},
	5752: {'tax_name' : 'Heterolobosea', 'rank' : 'phylum', 'count':0},
	33154: {'tax_name' : 'Opisthokonta', 'rank' : 'superkingdom_', 'count':0},
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
	"Unknown" : {'tax_name' : 'Unknown', 'rank' : 'NA', 'count':0},
	}, 
	'ranks' : {
	'superkingdom' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'superkingdom_' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'superkingdom__' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'kingdom' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'phylum' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'class' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'order' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'family' :{'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'genus' : {'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0},
	'species' :{'expected' : "NA", 'in_rank' : 0, 'not_in_rank' : 0}
	}}
	# now use taxa.csv to get the set of taxIDs nested under each 
	# lineage:
	
	return EntryDict

def init_test_bin_dict():
	# use the taxa.csv to make a dictionary with
	# all of the taxids in each of the tested lineages from entry dict
	TestBinDict = {554915 : {'tax_name' : 'Amoebozoa', 'rank' : 'superkingdom_', 'count':0},
		33154: {'tax_name' : 'Opisthokonta', 'rank' : 'superkingdom_', 'count':0},
		5878: {'tax_name' : 'Ciliophora', 'rank' : 'phylum', 'count':0},
		877183: {'tax_name' : 'Colpodellida', 'rank' : 'superkingdom___', 'count':0},
		3027: {'tax_name' : 'Cryptophyceae', 'rank' : 'class', 'count':0},
		2864: {'tax_name' : 'Dinophyceae', 'rank' : 'class', 'count':0},
		33682: {'tax_name' : 'Euglenozoa', 'rank' : 'phylum', 'count':0},
		#4751: {'tax_name' : 'Fungi', 'rank' : 'kingdom', 'count':0},
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
		10239: {'tax_name' : 'Viruses', 'rank' : 'superkingdom', 'count':0}}
	
	for test_id in TestBinDict.keys():
		test_rank = TestBinDict[test_id]['rank']
		bin_df = taxa_csv[taxa_csv[test_rank] == test_id]
		bin_ids = bin_df['tax_id'].tolist()
		TestBinDict[test_id]["nested_ids"] = bin_ids
	
	return TestBinDict

def get_expected_ranks(this_entry):
	# get the entry_id from this_entry:
	entry_id = int(this_entry.split('_')[0])
	# use metadata.csv to get the entry_taxid
	entry_taxid = metadat.loc[metadat['entry_id']==entry_id, 'tax_id'].item()
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
		EntryDict['tax_bin']["Unknown"]['count'] += 1
		return ""
	# check if the query_taxid is in TestBinDict
	for test_id in TestBinDict.keys():
		if query_taxid in TestBinDict[test_id]["nested_ids"]:
			EntryDict['tax_bin'][test_id]['count'] += 1
	# get the data for this tax_id from the taxa.csv file
	# taxon_dat = taxa_csv[taxa_csv['tax_id'] == query_taxid]
	# # iterate thru ranks:
	# for test_rank in test_ranks:
	# 	test_id = taxon_dat.loc[taxon_dat['tax_id']==query_taxid, test_rank].item()
	# 	try: test_id = int(test_id)
	# 	except ValueError: test_id = "NA"
	# 	if test_id in EntryDict['tax_bin'].keys():
	# 		# increment counter:
	# 		EntryDict['tax_bin'][test_id]['count'] += 1

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
	EntryDict['n_riboprotein'] = len(entry_data['aa_id'])
	EntryDict['n_core_pfams'] = len(core_pfam_list)
	expected_bin = EntryDict['expected_id']
	ingroup_count = "NA"
	outgroup_count = "NA"
	contam_pct = "NA"
	completeness = round(100*len(core_pfam_list) / len(ctg_pfam_ids63),2)
	if expected_bin != "None":
		ingroup_count = EntryDict['tax_bin'][expected_bin]['count'] 
		outgroup_count = 0
		for tax_bin in EntryDict['tax_bin'].keys():
			if tax_bin != expected_bin:
				outgroup_count += EntryDict['tax_bin'][tax_bin]['count']
		
		try: contam_pct = round((outgroup_count/(outgroup_count+ingroup_count)*100),2)
		except ZeroDivisionError: contam_pct = "NA"
	EntryDict['ingroup_count'] = str(ingroup_count)
	EntryDict['outgroup_count']	= str(outgroup_count)
	EntryDict['contam_pct'] = str(contam_pct)
	EntryDict['completeness'] = str(completeness) 

def init_output_file():
	# open the output file and write header
	outcsv_file = open(outcsv_path,'a')
	# output columns:
	out_cols1 = ["entry_handle", "entry_id", "tax_id", "n_seqs", "n_pfams", "tax_group", "contam_pct", "Amoebozoa", "Ciliophora", "Colpodellida", "Cryptophyceae", "Dinophyceae", "Euglenozoa", "Glaucocystophyceae", "Haptophyta", "Heterolobosea", "Opisthokonta", "Palpitomonas", "Perkinsozoa", "Rhizaria", "Rhodophyta", "Stramenopiles", "Viridiplantae", "Bacteria", "Archaea", "Viruses", "Other", "Unknown"]
	#out_cols2 = ['in_superkingdom', 'in_superkingdom_', 'in_superkingdom__', 'in_kingdom', 'in_phylum', 'in_class', 'in_order', 'in_family', 'in_genus', 'in_species']
	out_header = ",".join(out_cols1) + "\n"
	outcsv_file.write(out_header)
	return outcsv_file

def write_results(this_entry, EntryDict, outcsv_file):
	# write out the results for each entry
	out_line = [this_entry]
	out_line.append(str(EntryDict['entry_id']))
	out_line.append(str(EntryDict['tax_id']))
	out_line.append(str(EntryDict['n_riboprotein']))
	out_line.append(str(EntryDict['n_core_pfams']))
	out_line.append(str(EntryDict['expected_name']))
	out_line.append(str(EntryDict['contam_pct']))
	#out_line.append(str(EntryDict['completeness']))
	
	for taxbin in EntryDict['tax_bin'].keys():
		out_line.append(str(EntryDict['tax_bin'][taxbin]['count']))
	#for rank in EntryDict['ranks'].keys():
	#	out_line.append(str(EntryDict['ranks'][rank]['in_rank']))
	# write out:
	outline_ready = ",".join(out_line) + "\n"
	outcsv_file.write(outline_ready)

# input list of entry handles:
with open("mft_ribo_sequences.handles.txt") as f:
	entry_handles = f.read().splitlines()
	

outcsv_file = init_output_file()

# dict of ids to test for membership:
TestBinDict =  init_test_bin_dict()


# for each entry:
for this_entry in entry_handles:	
	print(this_entry)
	# skip if one of 28 removed entries:
	if this_entry in removed_entries:
		continue
	# gather and process the input files
	entry_data = collect_input_data(this_entry)
	# init entry dict:
	EntryDict = init_entry_dict(this_entry)
	# gather expected ranks:
	get_expected_ranks(this_entry)
	# iterate through pfams:
	entry_pfam_list = entry_data['pfam_id'].unique()
	# filter down list to set of core RPs:
	core_pfam_list = list(set(entry_pfam_list) & set(ctg_pfam_ids63))
	for this_pfam in core_pfam_list:
		process_entry_pfam(this_pfam, entry_data)
	# run statistics on expected and counted taxids
	entry_stats()
	# then write out EntryDict results:
	write_results(this_entry, EntryDict, outcsv_file)

# need this to write out output
outcsv_file.close()
