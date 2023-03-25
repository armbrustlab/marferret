import pandas as pd
import numpy as np

# Read in g1pa_tax_funct_euk_s, pfam_lookup, and TAXA_CSV.lin dataframes
g1pa_tax_funct_euk_s = pd.read_csv("g1pa.euk_s.tax_funct.csv")
pfam_lookup = pd.read_csv("Pfam-A.34.pfam_lookup.csv")
TAXA_CSV_lin = pd.read_csv("TAXA_CSV.lin.csv")

# Attach pfam_id and use that instead of name
g1pa_tax_funct_euk_s = pd.merge(g1pa_tax_funct_euk_s, pfam_lookup, on="pfam_name", how="left")

# Bring back tax_name, etc and recompute lineage
g1pa_tax_funct_euk_s = pd.merge(g1pa_tax_funct_euk_s, TAXA_CSV_lin[['tax_id', 'lineage']], on="tax_id", how="left")

# Merge stn and size to get a unique sample id
g1pa_tax_funct_euk_s['sample'] = g1pa_tax_funct_euk_s['stn'].astype(str) + '_' + g1pa_tax_funct_euk_s['size'].astype(str)

# Get unique tax_id-sample pairs and count number of stations
g1pa_tax_funct_euk_s_stn_tax = g1pa_tax_funct_euk_s[['tax_id', 'sample']].drop_duplicates().groupby('tax_id').agg(n_stn=('sample', 'nunique')).reset_index().sort_values(by='n_stn', ascending=False)

# Get tax_ids that are present in all 18 stations
all_stn_taxa = g1pa_tax_funct_euk_s_stn_tax[g1pa_tax_funct_euk_s_stn_tax['n_stn'] == 18]['tax_id'].unique()

# Filter g1pa_tax_funct_euk_s to only include taxa present in all 18 stations
g1pa_tax_funct_euk_s = g1pa_tax_funct_euk_s[g1pa_tax_funct_euk_s['tax_id'].isin(all_stn_taxa)]

# Display dimensions of g1pa_tax_funct_euk_s
print(g1pa_tax_funct_euk_s.shape)

# Get species * sample uniques and group to see frequency of taxa
g1pa_tax_funct_euk_s_stn_tax = g1pa_tax_funct_euk_s[['tax_id', 'sample']].drop_duplicates().groupby('tax_id').size().reset_index(name='n_stn').sort_values('n_stn', ascending=False)
print(g1pa_tax_funct_euk_s_stn_tax.shape)

# Keep only the taxa present in all 18 stations
g1pa_tax_funct_euk_s = g1pa_tax_funct_euk_s

# Merge pfam_id, lineage and create a unique sample id
g1pa_tax_funct_euk_s = pd.merge(g1pa_tax_funct_euk_s, pfam_lookup, on="pfam_name")
g1pa_tax_funct_euk_s = pd.merge(g1pa_tax_funct_euk_s, TAXA_CSV[["tax_id", "lineage"]], on="tax_id")
g1pa_tax_funct_euk_s["sample"] = g1pa_tax_funct_euk_s["stn"] + "_" + g1pa_tax_funct_euk_s["size"]

# Get tax_ids present in all 18 stations
all_stn_taxa = g1pa_tax_funct_euk_s.groupby('tax_id')['sample'].nunique()
all_stn_taxa = all_stn_taxa[all_stn_taxa == 18].index

# Filter g1pa_tax_funct_euk_s to only include taxa present in all 18 stations
g1pa_tax_funct_euk_s = g1pa_tax_funct_euk_s[g1pa_tax_funct_euk_s['tax_id'].isin(all_stn_taxa)]

# Create dataframe for core_pct with tax_id, n_pfams, and sample columns
g1pa_tax_funct_euk_s_core_pct = g1pa_tax_funct_euk_s[['tax_id', 'pfam_id']].drop_duplicates().groupby('tax_id').size().reset_index(name='n_pfams')
unique_samples = g1pa_tax_funct_euk_s['sample'].unique()
for sample in unique_samples:
    g1pa_tax_funct_euk_s_core_pct[sample] = None

# Function to compute pct_core for each species and sample
def compute_pct_core(species_df, ctg_dat, total_core):
    unique_samples = species_df['sample'].unique()
    pct_core_dict = {}
    for sample in unique_samples:
        species_sample_df = species_df[species_df['sample'] == sample][['tax_id', 'pfam_id']].drop_duplicates()
        n_core_sample = species_sample_df['pfam_id'].count()
        pct_core_dict[sample] = n_core_sample / total_core
    return pct_core_dict

# Iterate through species to compute and store pct_core
for species in g1pa_tax_funct_euk_s['tax_id'].unique():
    species_df = g1pa_tax_funct_euk_s[g1pa_tax_funct_euk_s['tax_id'] == species]
    this_lineage = species_df['lineage'].iloc[0]
    species_df = species_df[species_df['pfam_id'].isin(ctg_dat[ctg_dat['lineage'] == this_lineage]['pfam_id'].unique())]
    total_core = ctg_dat[ctg_dat['lineage'] == this_lineage]['pfam_id'].count()

    pct_core_dict = compute_pct_core(species_df, ctg_dat, total_core)
    for sample, pct_core in pct_core_dict.items():
        g1pa_tax_funct_euk_s_core_pct.loc[g1pa_tax_funct_euk_s_core_pct['tax_id'] == species, sample] = pct_core

# Check the dimension of the dataframe and display the first few rows
print(g1pa_tax_funct_euk_s_core_pct.shape)
print(g1pa_tax_funct_euk_s_core_pct.head())

# Save the results to a CSV file
g1pa_tax_funct_euk_s_core_pct.to_csv("data/g1pa/ctg/g1pa.euk_s.core_pct.per_sample.csv", index=False)
