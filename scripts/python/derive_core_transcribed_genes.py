import pandas as pd
from pathlib import Path
import argparse
import os


def handle_arguments():
	description = '''This script
		Example usage: ./combine_pfam_annotations.py -f taxIDs.txt MarFERReT.v1.Pfam34.annotations.csv
		'''
	parser = argparse.ArgumentParser(description=description)
	parser.add_argument("-m", "--metadata", help="metadata.csv file", type=str)
	parser.add_argument("-t", "--taxa_csv", help="taxa.csv file with NCBI Taxonomy relationships", type=str)
	parser.add_argument("-p", "--pfam_dat", help="Pfam annotation results in CSV file", type=str)
	parser.add_argument('output', type=str, help='Output CSV file with CTG results')

	return parser

def build_lineage_dict(lineages):
	# create a dictionary with taxonomy information for listed lineages
	LineagesDict = {}
	for this_lineage in lineages:
		lin_rank = taxa_csv[taxa_csv['tax_name'] == this_lineage]['rank'].item()
		lin_taxid = taxa_csv[taxa_csv['tax_name'] == this_lineage]['tax_id'].item()
		lin_taxa = taxa_csv[taxa_csv[lin_rank] == lin_taxid]['tax_id'].unique()
		LineagesDict[this_lineage] = {
		'tax_id' : lin_taxid, 
		'lineage_taxa' : lin_taxa,
		'rank' : lin_rank}
	return LineagesDict

def remove_nested_lineages(LineagesDict):
	# removes nested lineages
	for tax_name1 in LineagesDict.keys():
		tax_id1 = LineagesDict[tax_name1]['tax_id']
		lin_taxa1 = set(LineagesDict[tax_name1]['lineage_taxa'])	
		for tax_name2 in LineagesDict.keys():
			if tax_name1 != tax_name2 and tax_name2 != "Eukaryota":  
				tax_id2 = LineagesDict[tax_name2]['tax_id']
				lin_taxa2 = set(LineagesDict[tax_name2]['lineage_taxa'])			
				if tax_id1 in lin_taxa2:
					print("Overlap between ",tax_name1," & ",tax_name2)
					common_values = lin_taxa1.intersection(lin_taxa2)
					LineagesDict[tax_name2]['lineage_taxa'] = list(lin_taxa2 - common_values)
					print("Removed ",len(common_values)," taxa from ",tax_name2)
	return LineagesDict

def process_lineage(lineage, LineagesDict, pfam_pivot):
	# Filter data for the given lineage
	lineage_taxa = list(LineagesDict[lineage]['lineage_taxa'])
	#
	lineage_set = set(lineage_taxa) & set(pfam_pivot.columns)
	n_taxa = len(lineage_set)
	#
	pfam_sub = pfam_pivot[list(lineage_set)]
	#
	# Calculate rowsum and frequency
	pfam_sub['rowsum'] = pfam_sub.sum(axis=1)
	pfam_sub['frequency'] = pfam_sub['rowsum'] / n_taxa
	#
	# remove rows with frequency > 0:
	pfam_sub = pfam_sub[pfam_sub['frequency'] > 0]
	#
	# Filter core Pfam domains based on frequency
	core_pfams_subset = pfam_sub[pfam_sub['frequency'] >= 0.95].index
	#
	# Add lineage and pfam_id columns
	pfam_sub['pfam_id'] = pfam_sub.index
	pfam_sub['lineage'] = lineage
	ctg_dat_sub = pfam_sub[['lineage', 'rowsum', 'pfam_id', 'frequency']].rename(columns={'rowsum': 'n_taxa'})
	ctg_dat_sub = ctg_dat_sub[ctg_dat_sub['frequency'] >= 0.95]
	#
	# Save the results to a CSV file
	out_csv = mft_dir / f"ctg/MarFERReT.pfam_presence.{lineage}.ctg_catalog.v1.csv"
	pfam_sub.to_csv(out_csv, index=False)
	#
	return ctg_dat_sub 



def main():

	# parse arguments
	parser = handle_arguments()
	args = parser.parse_args()

	# Define directory paths
	mft_dir = Path("") # MarFERReT path
	data_dir = mft_dir / "data"
	pfam_dir = data_dir / "pfam"

	# Read input data
	mft_dat = pd.read_csv(data_dir / "MarFERReT.v1.metadata.csv")
	mft_dat = mft_dat.query('accepted == "Y"')
	mft_dat_trans = mft_dat.query('data_type in ["TSA", "SAT"]')
	taxa_csv = pd.read_csv(data_dir / "MarFERReT.v1.taxa.csv")
	pfam_dat = pd.read_csv(pfam_dir / "MarFERReT.v1.entry_pfam_sums.csv")

	# get set of unique Pfam domains
	unique_pfams = pfam_dat['pfam_name'].nunique()

	# Merge and filter data
	pfam_dat2 = (
		pfam_dat
		.merge(mft_dat[['entry_id', 'marferret_name', 'tax_id', 'data_type']], on='entry_id', how='left')
		.query('data_type in ["TSA", "SAT"]')
		#.merge(taxa_csv[['tax_id', 'species']], on='tax_id', how='left')
		.drop(['entry_id', 'marferret_name', 'tax_id', 'data_type'], axis=1)
	)

	# presence/absence table:
	pfam_dat_matched = pfam_dat2[['pfam_id','tax_id']]
	pfam_dat_matched = pfam_dat_matched.drop_duplicates()
	pfam_dat_matched.shape
	pfam_dat_matched['observed'] = 1

	pfam_pivot = (
	  pfam_dat_matched
	  .pivot_table(index='pfam_id', columns='tax_id', values='observed', fill_value=0))

	# List the lineages to define core genes for
	lineages = [
		'Dinophyceae',
		'Bacillariophyta',
		'Ciliophora',
		'Haptophyta',
		'Chlorophyta',
		'Rhizaria',
		'Ochrophyta',
		'Amoebozoa',
		'Opisthokonta',
		'Eukaryota'
	]

	# Create the lineage dictionary
	LineagesDict = build_lineage_dict(lineages)

	# Remove nested lineages if desired;
	# checks if lineages overlap and removes ID sets except for Eukaryota
	LineagesDict = remove_nested_lineages(LineagesDict)

	for lineage in LineagesDict.keys():
		print(lineage," : ",len(LineagesDict[lineage]['lineage_taxa']))

	# Process each lineage and save the results
	ctg_dat = pd.concat([process_lineage(lineage, LineagesDict, pfam_pivot) for lineage in lineages])
	ctg_dat.to_csv(mft_dir / "ctg/MarFERReT.core_genes.v1.csv", index=False)

if __name__ == "__main__":
	main()
