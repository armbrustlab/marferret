# Dr. Ryan D. Groussman
# Armbrust Lab
# University of Washington, 2023

# Zenodo repository for MARMICRODB:
# https://zenodo.org/record/3520509
# MARMICRODB database for taxonomic classification of (marine) metagenomes

# These files were downloaded locally from Zenodo on 04/15/21
wget https://zenodo.org/record/3520509/files/MARMICRODB.faa.bz2
wget https://zenodo.org/record/3520509/files/MARMICRODB.fmi
wget https://zenodo.org/record/3520509/files/MARMICRODB_catalog.tsv
wget https://zenodo.org/record/3520509/files/names.dmp
wget https://zenodo.org/record/3520509/files/nodes.dmp
wget https://zenodo.org/record/3520509/files/phylogenies.tar.gz
wget https://zenodo.org/record/3520509/files/scripts.tar.gz

# unpack:
bzip2 -dk MARMICRODB.faa.bz2

# Process the MARMICRODB fasta file to keep an entry only if:
	# 1. it is not a MAG
	# 2. it is not a eukaryote,
	# 3. it has a valid NCBI taxID

# MarMicroDB fasta headers have the defline format:
	# >($1)-($2)_$3
# where,
	# $1 = 'genome' id
	# $2 = seq_n
	# $3 = tax_id

# The 'MARMICRODB_catalog.tsv' has metadata about the data in MARMICRODB.
# The following code uses this file to identify the MARMICRODB entries
# that satisfy our three criteria.

## Initiate an R session ##

# import dependencies
library(tidyverse)

# You should set this to wherever you downloaded the MARMICRODB data.
setwd("/mnt/nfs/projects/ryan/EukRefDB_2021/MARMICRODB")

# import the metadata file
mmdb_cat = read.table("MARMICRODB_catalog.tsv", sep="\t", quote = "\"", header = TRUE, stringsAsFactors = FALSE)
# there are 18768 entries in MARMICRODB:
dim(mmdb_cat) # "[1] 18768    14"

# we want to filter first by sequence type:
mmdb_cat$sequence_type %>% unique()
	# [1] "mag"           "isolate"       "sag"           "transcriptome"

# start a new df that exclused the metagenome assembled genomes (mag):
mmdb_cat_kept = mmdb_cat %>% filter(sequence_type != "mag")
# this trims the entry size down by more than half, but
# now satisfies a major part of our data product needs.
dim(mmdb_cat_kept) # "[1] 8106   14"

# now we filter it by taxonomic domain.
mmdb_cat_kept$domain %>% unique()
	# [1] "archaea"   "bacteria"  "eukaryote"

# The eukaryotes in MARMICRODB are largely MMETSP species and are
# redundant to the ones in MarFERReT, so we can remove them:
mmdb_cat_kept = mmdb_cat_kept %>% filter(domain != "eukaryote")
# this only removes 185 entries; a small portion of the library.
dim(mmdb_cat_kept) # [1] 7921   14

# now we filter out entries that do not have an NCBI taxID:
mmdb_cat_kept = mmdb_cat_kept %>% filter(!is.na(taxid))
# this is a small reduction; only 17 entries are removed.
dim(mmdb_cat_kept)
"[1] 7904   14"

# There were 53 additional entries have NCBI taxIDs that
# were not compatible with the NCBI Taxonomy database at the time of analysis.
# There were also removed.
remove_ids=c('1046724', '1054037', '1104325', '1109743', '1118061', '1121105', '1122931', '1137271', '1156433', '1217713', '1220582', '1223307', '1232666', '1267535', '1313292', '1317118', '1331660', '140100', '1417988', '1449126', '1453999', '1454004', '1457154', '1609981', '1915400', '319236', '335659', '338187', '440512', '443255', '469595', '469596', '595593', '62153', '62928', '640511', '644968', '656024', '667632', '67281', '68178', '68201', '684949', '741091', '861208', '861530', '876044', '911239', '913102', '935543', '941770', '97139', '981223')
length(remove_ids) # [1] 53
mmdb_cat_kept = mmdb_cat_kept %>% filter(!taxid %in% remove_ids)
# dim(mmdb_cat_kept) "[1] 7851   14"

# the 'genome' field in the metadata csv is the unique identifier
# for each entry; so we'll output a list of these and the taxID
# for use in a python script to extract sequences in these entries
# from the main fasta file.
mmdb_cat_kept %>% pull(genome) %>% write.table("/mnt/nfs/projects/ryan/EukRefDB_2021/MARMICRODB/MARMICRODB_catalog.filtered.genome.txt", quote = FALSE, col.names = FALSE, row.names	= FALSE)

# We're done with R now
quit()

# Outside of R, use bash to execute a python script.
# Run a custom python script to keep the sequences for the entry IDs found in
# the 'MARMICRODB_catalog.filtered.genome.txt'

# Base directories:
MMDB_DIR="/mnt/nfs/projects/ryan/EukRefDB_2021/MARMICRODB"
MARFERRET_DIR="/scratch/ryan/mft/marferret-main/"

# path to the big MARMICRODB file:
MARMICRODB_FASTA="${MMDB_DIR}/MARMICRODB.faa"
KEPT_IDS="${MMDB_DIR}/MARMICRODB_catalog.filtered.genome.txt"

# The output is placed in the MarFERReT directory for incorporation with diamond:
OUT_FASTA="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered.uid.faa"
MARMICRODB_UID2TAX="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered.uid2tax.tab"
EukRefDB.trim_marmicrodb.py -k ${KEPT_IDS} -f ${MARMICRODB_FASTA}  -o ${OUT_FASTA} -t ${MARMICRODB_UID2TAX}
# This outputs a filtered fasta file (OUT_FASTA)
# and a uid2tax file for DIAMOND annotation (MARMICRODB_UID2TAX)

# For reasons unknown to the author, some MARMICRODB sequences
# have numeric values in the amino acid sequence field:
# for example:
"	>mmdb1 1A-1_54252
	21826632182675MHGKSGSRGGPVAQPGRALGSHPRGPGFKSRPVHHPLLDALRAWGEL
	>mmdb2 1A-2_54252
	21826632182675MDVLDEVFERVVKARIFRNRSVLSPDYIPDKLPHREREIRALGSIV"

# This script removes the sequences with numeric values:
MARMICRODB_FASTA="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered.uid.faa"
OUTPUT_FASTA="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered2.uid.faa"
MarMicroDB.remove_bad_seqs.py -f ${MARMICRODB_FASTA} -o ${OUTPUT_FASTA}

# compare:
grep -c ">" $MARMICRODB_FASTA # 27890788
grep -c ">" $OUTPUT_FASTA # 27566953
# This removes 323835 sequences with numerics (1.2%)

# compress the MMDB uid2tax and FASTA file:
gzip $MARMICRODB_FASTA
gzip $MARMICRODB_UID2TAX

### Merging with MarFERReT

# concatenate the merged files
cd ${MARFERRET_DIR}/data/marmicrodb/
MARFERRET_SEQS="${MARFERRET_DIR}/data/MarFERReT.v1.proteins.faa.gz"
MARMICRODB_SEQS="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered2.uid.faa.gz"

ls -lht $MARFERRET_SEQS # 4.8G
ls -lht $MARMICRODB_SEQS # 5.0G

cat $MARFERRET_SEQS $MARMICRODB_SEQS > MarFERReT.v1.MMDB.combined.faa.gz
COMBINED_SEQS="${MARFERRET_DIR}/data/marmicrodb/MarFERReT.v1.MMDB.combined.faa.gz"
ls -lht $COMBINED_SEQS # 9.8G
zgrep -c ">" $COMBINED_SEQS # 58368183

# combine both uid2tax files
# the order matters here; the MarFERReT.v1.taxonomies.tab.gz contains
# the necessary header, and the marmicrodb.filtered.uid2tax.tab.gz
# file does not contain a header so that it can be added right on.
MARFERRET_UID2TAX="${MARFERRET_DIR}/data/MarFERReT.v1.taxonomies.tab.gz"
MARMICRODB_UID2TAX="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered.uid2tax.tab.gz"
cat $MARFERRET_UID2TAX $MARMICRODB_UID2TAX > MarFERReT.v1.MMDB.combined.uid2tax.tab.gz
ls -lht MarFERReT.v1.MMDB.combined.uid2tax.tab.gz # 165M

cd ${MARFERRET_DIR}/data/marmicrodb/dmnd/

#### Build the DIAMOND database
NCORES=16 # adjust this value for your system

COMBINED_SEQS="${MARFERRET_DIR}/data/marmicrodb/MarFERReT.v1.MMDB.combined.faa.gz"
TAXONNODES="${MARFERRET_DIR}/data/diamond/ncbi/nodes.dmp"
TAXONNAMES="${MARFERRET_DIR}/data/diamond/ncbi/names.dmp"
TAXONMAP="${MARFERRET_DIR}/data/marmicrodb/MarFERReT.v1.MMDB.combined.uid2tax.tab.gz"
MFT_MMDB_DMND_DB="${MARFERRET_DIR}/data/marmicrodb/dmnd/MarFERReT.v1.MMDB.combined.dmnd"
time diamond makedb --threads ${NCORES} --in $COMBINED_SEQS --db ${MFT_MMDB_DMND_DB} --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} --taxonmap ${TAXONMAP}
