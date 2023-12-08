Every MarFERReT receives an [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) identifier (taxID) through inheritance of a taxID from the source metadata or manual assignment of a taxID, if not existing assignment was available. All entries were manually updated if necessary to reflect changes and additions to the NCBI Taxonomy database, as of October 11th, 2022.

NCBI Taxonomy in whole is a large database, so we like to use the [‘taxtastic’ v0.9.2 software package](https://github.com/fhcrc/taxtastic) to generate a smaller table (we reference these as 'taxa.csv' files) that contains the name, rank, and nested taxonomic classifications for a set of given NCBI taxIDs. [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) taxIDs are available in the `metadata.csv` file available in the Zenodo data repository. 

Here is how we created a `taxa.csv` file from the command line:
``` shell

# Define directory structures:
MARFERRET_DIR=""

# extract the taxIDs from MarFERReT.v1.metadata.csv as a list:
TAXID_LIST="${MARFERRET_DIR}/data/MarFERReT.v1.taxIDs.txt" 
# (don't include the header from the metadata.csv file)

# Ensure the NCBI information is downloaded here
cd ${MARFERRET_DIR}/data/diamond/ncbi

# now use the 'taxtable' function with 'taxit' to build the taxa.csv file:
taxit taxtable -f MarFERReT.v1.taxIDs.txt -o MarFERReT.v1.taxa.csv ${MARFERRET_DIR}/data/diamond/ncbi/ncbi_taxonomy.db

# This outputs the `MarFERReT.v1.taxa.csv` file:
MARFERRET_TAXA_CSV="MarFERReT.v1.taxa.csv"
```

Different `taxa.csv` files can be built for different combinations of taxIDs. We generate taxIDs for the combined MarFERReT and MARMICRODB database detailed here

We generate a large `taxa.csv` file for all of the taxIDs represented in the downloaded sequences associated with ribosomal protein Pfam families for the RP63 analysis here:

Here is how we generated the combined MarFERReT + MARMICRODB taxa.csv:
``` shell

#### MarFERReT_v1_MMDB TAXA_CSV ####

# MarFERReT taxIDs:
MARFERRET_IDS="${MARFERRET_DIR}/data/MarFERReT.v1.taxIDs.txt"
NCBI_DIR="${MARFERRET_DIR}/data/diamond/ncbi"
# path to collected MARMICRODB IDs:
MMDB_IDS="${MARFERRET_DIR}/data/marmicrodb/MARMICRODB_catalog.filtered.taxid.txt"

cd ${MARFERRET_DIR}/data/marmicrodb
# combine the unique taxIDs for MarFERReT and MarMicroDB filtered:
cat ${MARFERRET_IDS} ${MMDB_IDS} | sort | uniq > MarFERReT_v1_MARMICRODB.taxids.txt

# build:
taxit taxtable -f MarFERReT_v1_MARMICRODB.taxids.txt -o MarFERReT_v1_MMDB.taxa.csv ${NCBI_DIR}/ncbi_taxonomy.db

# this produces the file: MarFERReT_v1_MMDB.taxa.csv

```
