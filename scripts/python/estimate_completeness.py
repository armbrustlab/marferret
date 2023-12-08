import os
import pandas as pd

# Relative path for MarFERReT installation:
mft_dir = ""
os.chdir(mft_dir)

def assign_lineage(row):
    if row['class'] == '2864':
        return 'Dinophyceae'
    elif row['phylum'] == '2836':
        return 'Bacillariophyta'
    elif row['phylum'] == '5878':
        return 'Ciliophora'
    elif row['phylum'] == '2830':
        return 'Haptophyta'
    elif row['phylum'] == '3041':
        return 'Chlorophyta'
    elif row['superkingdom__'] == '543769':
        return 'Rhizaria'
    elif row['superkingdom___'] == '2696291':
        return 'Ochrophyta'
    elif row['superkingdom_'] == '554915':
        return 'Amoebozoa'
    elif row['superkingdom_'] == '33154':
        return 'Opisthokonta'
    elif row['superkingdom'] == '2759':
        return 'Eukaryota'
    else:
        return 'Other'

# Load data
ctg_dat = pd.read_csv(os.path.join(mft_dir, "data/ctg/MarFERReT.core_genes.v1.csv"), dtype=str)
g1pa_tax_funct_sums = pd.read_csv(os.path.join(mft_dir, "data/g1pa/compare_db.mft_mmdb_pfam_n.csv"), dtype=str)
TAXA_CSV = pd.read_csv(os.path.join(mft_dir, "data/MarFERReT_v1_MMDB.taxa.csv"), dtype=str)

# Assign lineage to taxa
TAXA_CSV['lineage'] = TAXA_CSV.apply(assign_lineage, axis=1)

# Merge PFAM data
g1pa_tax_funct_sums = g1pa_tax_funct_sums.merge(TAXA_CSV[['tax_id', 'rank', 'tax_name', 'species']], on='tax_id', how='left')

# Group by pfam and species, and calculate sums
species_taxfun = g1pa_tax_funct_sums[g1pa_tax_funct_sums['species'].notna()].groupby(['species', 'pfam_name']).agg({'n': 'sum'}).reset_index()
species_taxfun.rename(columns={'species': 'tax_id'}, inplace=True)

# Merge taxonomic data
species_taxfun = pd.merge(species_taxfun, TAXA_CSV[['tax_id', 'tax_name', 'superkingdom']], on='tax_id', how='left')

# Restrict to eukaryotes only
species_taxfun = species_taxfun[species_taxfun['superkingdom'] == '2759']

# Attach pfam_id
pfam_lookup = pd.read_csv("Pfam-A.34.pfam_lookup.csv")
species_taxfun = pd.merge(species_taxfun, pfam_lookup, on='pfam_name', how='left')

# Add lineage to species_taxfun
species_taxfun = pd.merge(species_taxfun, TAXA_CSV[['tax_id', 'class', 'phylum', 'superkingdom_', 'superkingdom__', 'superkingdom___']], on='tax_id', how='left')
species_taxfun['lineage'] = species_taxfun.apply(assign_lineage, axis=1)

# Drop unnecessary columns
species_taxfun.drop(columns=['class', 'phylum', 'superkingdom_', 'superkingdom__', 'superkingdom___'], inplace=True)

# Group by tax_id to get total contigs, tax_name, lineage, and number of unique pfams per species
species_sums = species_taxfun.groupby('tax_id').agg(
    n_contigs=('n', 'sum'),
    tax_name=('tax_name', 'first'),
    lineage=('lineage', 'first'),
    n_pfams=('pfam_name', 'nunique')
).reset_index().sort_values(by='n_contigs', ascending=False)

# Add columns for n_core and total_core, set default values to NA
species_sums['n_core'] = pd.NA
species_sums['total_core'] = pd.NA

# Function to calculate completeness
def calculate_completeness(this_species, species_taxfun, ctg_dat):
    species_df = species_taxfun[species_taxfun['tax_id'] == this_species]
    this_lineage = species_df['lineage'].iloc[0]
    species_df = species_df[species_df['pfam_id'].isin(ctg_dat[ctg_dat['lineage'] == this_lineage]['pfam_id'].unique())]
    n_core = len(species_df['pfam_id'])
    total_core = len(ctg_dat[ctg_dat['lineage'] == this_lineage]['pfam_id'].unique())
    return n_core, total_core

# Iterate through each species to calculate completeness
for this_species in species_taxfun['tax_id'].unique():
    n_core, total_core = calculate_completeness(this_species, species_taxfun, ctg_dat)
    species_sums.loc[species_sums['tax_id'] == this_species, 'n_core'] = n_core
    species_sums.loc[species_sums['tax_id'] == this_species, 'total_core'] = total_core

# Calculate pct_core
species_sums['pct_core'] = species_sums['n_core'] / species_sums['total_core']

# Display completeness statistics
print("Max pct_core: ", species_sums['pct_core'].max())
print("Min pct_core: ", species_sums['pct_core'].min())
print("Mean pct_core: ", species_sums['pct_core'].mean())

# Save species_sums as a csv file
species_sums.to_csv("MarFERReT_v1.g1pa.pct_core.csv", index=False)
