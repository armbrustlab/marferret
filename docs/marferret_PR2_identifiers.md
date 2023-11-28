MarFERReT entries have been linked to matching entries in the PR2 database ecosystem of protist ribosomal reference sequences (https://pr2-database.org/), along with 9 classification levels from the community curated PR2 database. We linked the MarFERReT entries to PR2 entries using shared NCBI taxID information and preserved strain-level diversity information using fuzzy string matching to identify the best matching PR2 taxa at the strain- or species- level wherever possible. 

This process can be broken down into three major steps:

1. Download the PR2 taxonomy database files and retrieve NCBI taxonomy information using Entrez API

Access the PR2 database release from source from here:
https://github.com/pr2database/pr2database/releases

We use PR2 version 5.0.0, and accessed the full release table (pr2_version_5.0.0_merged.xlsx) on August 16th, 2023. 

We use a custom python script to python to retrieve the NCBI taxID information linked to the GenBank accession identifier for every PR2 sequence 
``` python 
import pandas as pd
import time
import json
from Bio import Entrez
import os 

# add your email address here:
Entrez.email = '' 
# obtain an API key to speed up the retrieval process
Entrez.api_key = '' 

def accession2taxid(acc: str, db="nucleotide") -> str:
    handle = Entrez.esearch(db=db, term=acc)
    record = Entrez.read(handle)
    gi = record["IdList"][0]
    handle = Entrez.esummary(db=db, id=gi, retmode="json")
    result = json.load(handle)["result"]
    taxid = result[gi]["taxid"]
    title = result[gi]["title"].replace(",","")
    return str(taxid), str(title)

# Load the full PR2 table 
pr2_dat = pd.read_excel("pr2_version_5.0.0_merged.xlsx")

# we want to retain these columns:
subset_columns = ['pr2_accession', 'genbank_accession', 'domain', 'supergroup', 'division', 'subdivision', 'class', 'order', 'family', 'genus', 'species'] 
pr2subset = pr2_dat[subset_columns].copy()

# we only need to use the pr2 accession and genbank accession for the next step:
pr2_list1 = pr2subset[['pr2_accession','genbank_accession']]

# the next step is to use the 'genbank_accession' to get the NCBI taxID
outdir = ""
outfilename = "pr2.ncbi_tax.csv" 
output_file_path = os.path.join(outdir, outfilename)

# Open the file for writing at the beginning 
with open(output_file_path, 'w') as f:
       f.write(",".join(pr2_list1.columns) + ",tax_id,title" + "\n")

iteration = 0

# Loop through the dataframe with genbank_accession
for index, row in pr2_list2.iterrows():
       print(row['genbank_accession'])
       gb = row['genbank_accession']
       pr2 = row['pr2_accession']
       try:
              taxid, title = accession2taxid(gb, db="nucleotide")  
       except Exception as e:
              print(f"Error: {e}")
              continue  # Continue to the next iteration
       with open(output_file_path, 'a') as f:
              f.write(",".join([pr2, gb, taxid, title]) + "\n")
       iteration += 1
       print(iteration)

# merge this data back with the original dataframe:
pr2_taxid_dat = pd.read_csv("pr2.ncbi_tax.csv")
subcols = ['pr2_accession', 'tax_id', 'title']

pr2_df_v5_ncbi = pd.merge(pr2_df_v5, pr2_taxid_dat[subcols], on='pr2_accession', how='left')
pr2_df_v5_ncbi.shape # (221085, 13)

# save out the PR2 data with NCBI taxIDs added:
pr2_df_v5_ncbi.to_csv("pr2_version_5.ncbi_taxid.csv", index=False)
```

2. Subset the PR2 database using the NCBI taxIDs from MarFERReT entries

Next, the pool of PR2 IDs and MarFERReT entries are linked using shared NCBI taxID information: 
``` python
mft_dat = pd.read_csv("data/MarFERReT.v1.metadata.csv")

# select cols to keep for the pr2merge file:
mft_cols = ['entry_id', 'tax_id', 'marferret_name', 'original_taxID', 'source_name', 'alias'] 

# merge with marferret on 'taxID'
mft_pr2_merge = pd.merge(mft_dat[mft_cols], pr2_df_v5_ncbi.dropna(subset=['tax_id'], on='tax_id', how='left')

# save out intermediate:
mft_pr2_merge.to_csv("MarFERReT.v1.mft_pr2_merge.csv.gz", index=False)
```


3. Use fuzzy string matching to suggest the best match among shared  IDs

 A custom python script was used with the updated NCBI taxIDs associated with MarFERReT entries to identify the sequence tags in the PR2 database that share the same NCBI taxID. We used RapidFuzz, a rapid fuzzy string matching algorithm to help identify the PR2 sequence tags matching MarFERReT entries at the strain-level wherever possible, as the strain identities are not always captured in the NCBI Taxonomy framework and the PR2 taxonomy levels only descend to the species level, even if strain information is included in the ID description. 
 
We wrote this custom python script to output the top PR2 matches to MarFERReT entries from the set of PR2 identifiers sharing the same NCBI taxID, here: [MarFERReT.NCBI_to_PR2.py ](https://github.com/armbrustlab/marferret/blob/main/scripts/python/MarFERReT.NCBI_to_PR2.py)

After the script above generates top-matching PR2 entries to MarFERReT entries, each entry was manually reviewed to identify the most accurate strain-level match, proceeding to species-level matches if no strain match was found, and genus-level if no species-level match exists in the PR2 database. These final selected PR2 identifiers are listed along with PR2 taxonomy information in MarFERReT.v1.metadata.csv.
