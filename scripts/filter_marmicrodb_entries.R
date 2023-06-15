#!/usr/bin/env R
# AUTHOR: Ryan Groussman

## Initiate an R session ##

# import dependencies
library(tidyverse)

# import the MARMICRODB metadata file
mmdb_cat = read.table("MARMICRODB_catalog.tsv", sep="\t", quote = "\"", header = TRUE, stringsAsFactors = FALSE)
# there are 18768 entries in MARMICRODB:
dim(mmdb_cat) # "[1] 18768    14"

# we want to filter first by sequence type:
mmdb_cat$sequence_type %>% unique()
# [1] "mag"           "isolate"       "sag"           "transcriptome"

# start a new df that excludes the metagenome assembled genomes (mag):
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
dim(mmdb_cat_kept) # "[1] 7904   14"

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
mmdb_cat_kept %>% pull(genome) %>% write.table("MARMICRODB_catalog.filtered.genome.txt", quote = FALSE, col.names = FALSE, row.names	= FALSE)

# Conclude the R session
quit()
