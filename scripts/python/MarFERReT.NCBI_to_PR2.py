#!/usr/bin/env python3

import os
import pandas as pd

# import marferret data for each entry:
# now lets bring in the MarFERReT (full csv table)
mft_dat = pd.read_csv("MarFERReT.v1.entry_curation.csv")

# select the marferret cols we want for the pr2merge file:
mft_cols = ['candidate_id', 'tax_id', 'marferret_name', 'original_taxID', 'source_name', 'alias'] 

# iterate thru the df and build a dict:
MftDict = dict()

for index, row in mft_dat.iterrows():
	id = row['candidate_id']
	MftDict[id] = dict()
	for col in ['marferret_name', 'source_name', 'alias']:
		MftDict[id][col] = row[col]
	MftDict[id]['tax_id'] = int(row['tax_id'])
	# concatenated strings:
	MftDict[id]['query_str'] = str(row['marferret_name']) + ' ' + str(row['source_name']) + ' ' + str(row['alias'])

len(MftDict) # 902

# import PR2 data:
PR2toNCBI_dat = pd.read_csv("pr2_version_5.ncbi_taxid.csv.gz")
PR2toNCBI_dat.shape # (221085, 13)

# drop PR2 entries without an NCBI taxID:
PR2toNCBI_dat = PR2toNCBI_dat.dropna(subset=['tax_id'])
PR2toNCBI_dat['tax_id'] = PR2toNCBI_dat['tax_id'].astype(int)
PR2toNCBI_dat.shape # (220871, 13)

# filter down to taxIDs in the MarFERReT records:
marferret_taxids = set(mft_dat['tax_id'])
len(marferret_taxids) # 516

# remove these two from the record:
"""candidate_id	tax_id	marferret_name	original_taxID	source_name	alias	uniq_pr2_ct	uniq_gbacc_ct
540	100272	unclassified_eukaryote	100272	Unidentified sp. Strain CCMP1999	MMETSP1475	28238	28000
616	187299	uncultured_ciliate	187299	Undescribed Undescribed Strain Undescribed	MMETSP1317	868	868"""
marferret_taxids.remove(100272)
marferret_taxids.remove(187299)
len(marferret_taxids) # 514

PR2toNCBI_dat = PR2toNCBI_dat[PR2toNCBI_dat['tax_id'].isin(marferret_taxids)]
PR2toNCBI_dat.shape # (4579, 13)

# Create the PR2toNCBI dictionary with these columns:
pr2_cols = ['tax_id', 'pr2_accession', 'genbank_accession', 'species', 'title']

# filter the PR2Dict to remove in 'species'
remove_species_with = [':plas']
# remove choices from consideration with these words in 'title':
remove_titles_with = ['16s', 'chloroplast', 'plastid']

PR2Dict = dict()

for index, row in PR2toNCBI_dat.iterrows():
	# remove species that contains: ":plas" 
	if ":plas" in row['species']:
		continue
	# skip titles with these substrings:
	if any(substring in row['title'] for substring in remove_titles_with):
		continue
	tax_id = row['tax_id']
	if tax_id not in PR2Dict.keys():
		PR2Dict[tax_id] = dict()
	pr2_accession = row['pr2_accession']
	PR2Dict[tax_id][pr2_accession] = dict()
	PR2Dict[tax_id][pr2_accession]['pr2_accession'] = pr2_accession
	PR2Dict[tax_id][pr2_accession]['species'] = row['species']
	PR2Dict[tax_id][pr2_accession]['title'] = row['title']
	# Concatenate the columns 'title' and 'species'
	PR2Dict[tax_id][pr2_accession]['choice_str'] = row['title'] + ' ' + row['species']


# from the PR2Dict
# 	1. make the choice string
# 	2. process title / accession to get 
# 	3. remove, filter, clean, etc
# 	4. get unique values in list (does it disrupt order?)

from rapidfuzz.process import extractOne
from rapidfuzz.distance import Levenshtein
from rapidfuzz.fuzz import token_set_ratio
from rapidfuzz.fuzz import token_sort_ratio
from rapidfuzz.fuzz import token_ratio
from rapidfuzz import fuzz, utils

# common tokens to remove in the 99% highest freq
remove_tokens = ['18s', 'and', 'clone', 'collection', 'complete', 'contig', 'cs', 'culture', 'dcm', 'eukaryote', 'for', 'gene', 'genome', 'genomic', 'internal', 'isolate', 'large', 'partial', 'ribosomal', 'rna', 'rrna', 'sequence', 'shotgun', 'small', 'spacer', 'strain', 'subunit', 'transcribed', 'type', 'uncultured', 'whole']

def token_ratio_by_taxID(query_str, pr2_subdict):
	# clean query string:
	query_str = utils.default_process(query_str)
	# iterate through pr2_subdict
	for pr2_id in pr2_subdict.keys():
		choice_str = pr2_subdict[pr2_id]['choice_str']
		# clean and process choice string:
		choice_str = utils.default_process(choice_str)
		# Remove abundant non-informative tokens:
		for token in remove_tokens:
			choice_str = choice_str.replace(token, '')
		# Determine similarity and assign to dictionary
		pr2_subdict[pr2_id]['similarity'] = fuzz.token_ratio(query_str,choice_str)
		
	# sort pr2_subdict based on the 'similarity' and retrieve top three
	sorted_dict = sorted(pr2_subdict.items(), key=lambda x: x[1]['similarity'], reverse=True)
	
	# return the top three PR2 entries; merging pr2_accession|species|title
	top3_pr2 = sorted_dict[:3]
	top3_list = []
	for n in range(len(top3_pr2)):
		top3_list.append("|".join([top3_pr2[n][1]['pr2_accession'],top3_pr2[n][1]['species'],top3_pr2[n][1]['title']]))
	return top3_list


# iterate through each entry of the marferret dict:
for mft_id in MftDict.keys():
	print("Entry #: ",mft_id)
	# tax_id for the marferret entry:
	tax_id = MftDict[mft_id]['tax_id']
	# set initial state for matches:
	MftDict[mft_id]['match1'] =  "None"
	MftDict[mft_id]['match2'] =  "None"
	MftDict[mft_id]['match3'] =  "None"
	# check if this tax_id is in PR2Dict
	if tax_id in PR2Dict.keys():
		# if true, run fuzzy match on PR2 entries with matching taxID
		# return top three matches (concatenated )
		query_str = MftDict[mft_id]['query_str']
		print("Query #: ",query_str)
		top3_matches = token_ratio_by_taxID(query_str, PR2Dict[tax_id])
		print("Top: ",top3_matches[0])
		MftDict[mft_id]['match1'] = top3_matches[0]
		if len(top3_matches) > 1:
			MftDict[mft_id]['match2'] = top3_matches[1]
		if len(top3_matches) > 2:
			MftDict[mft_id]['match3'] = top3_matches[2]
	else:
		# if not, run fuzzy match on family-level
		print("Query not in PR2 dict")

# now write out results and evaluate the matching:
out_file_path = "MarFERReT.v1.mft_pr2_matches.csv"
headers = ",".join(["candidate_id","tax_id","marferret_name","source_name","alias","match1","match2","match3"])

with open(out_file_path, 'a') as f:
	f.write(headers + "\n")
	# output columns: 
	for mft_id in sorted(MftDict.keys()):
		lineout = [str(mft_id)]
		for value in ["tax_id","marferret_name","source_name","alias","match1","match2","match3"]:
			lineout.append(str(MftDict[mft_id][value]))
		f.write(",".join(lineout)+"\n")
	


