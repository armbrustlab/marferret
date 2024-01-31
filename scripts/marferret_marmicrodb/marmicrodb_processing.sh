#!/bin/bash

# AUTHOR: Ryan Groussman, PhD

# PURPOSE: This script contains bash shell code to download and 
# perform pre-processing steps on the MARMICRODB v1.0 sequence
# database prior to merge with the MarFERReT eukaryotic library

# 1. Download MARMICRODB database from the Zenodo repository
# 2. Filter out MAGs, eukaryotes, and entries without NCBI taxIDs
# 3. Remove a small subset of sequences with numerical values in the sequence field
# 4. Run a script to generate a FASTA file and taxonomy table for combining with MarFERReT



# MARMICRODB reference:
# Hogle, S. L. MARMICRODB database for taxonomic classification of (marine) metagenomes (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.3520509 (2019). 

# Zenodo repository for MARMICRODB:
# https://zenodo.org/record/3520509

# MICRODB data products can be downloaded from Zenodo directly from the command line using wget:
wget https://zenodo.org/record/3520509/files/MARMICRODB.faa.bz2
wget https://zenodo.org/record/3520509/files/MARMICRODB.fmi
wget https://zenodo.org/record/3520509/files/MARMICRODB_catalog.tsv
wget https://zenodo.org/record/3520509/files/names.dmp
wget https://zenodo.org/record/3520509/files/nodes.dmp
wget https://zenodo.org/record/3520509/files/phylogenies.tar.gz
wget https://zenodo.org/record/3520509/files/scripts.tar.gz

# unzip the FASTA file:
bzip2 -dk MARMICRODB.faa.bz2

# For our goal of integrating MARMICRODB with the eukaryotic MarFERReT library,
# we are removing some of the entries from MARMICRODB:
	# 1. Metagenome-assembled genomes (MAGs), for stringent reference organism identity
	# 2. MMETSP marine eukaryotes, which are redundant with MarFERReT sequence data.
	# 3. Entries without a valid NCBI taxID

# For this, we want to process the MARMICRODB fasta file to keep an entry only if:
	# 1. it is not a MAG
	# 2. it is not a eukaryote,
	# 3. it has a valid NCBI taxID

# The FASTA defline in MARMICRODB.faa has:
	# >($1)-($2)_$3
# where,
	# $1 = 'genome' id
	# $2 = seq_n
	# $3 = tax_id

# The 'MARMICRODB_catalog.tsv' has metadata about the data in MARMICRODB.
# The following code uses this file to identify the MARMICRODB entries
# that satisfy our three criteria.

# Run this R script to filter out MAGs, eukaryotes, and entries without NCBI taxIDs:
# this script takes as input the 'MARMICRODB_catalog.tsv' file and outputs a file called
# 'MARMICRODB_catalog.filtered.genome.txt', containing the filtered 'genome' IDs from MARMICRODB.
./filter_marmicrodb_entries.R

# Outside of R, use bash to execute a python script.
# Run a custom python script to keep the sequences for the entry IDs found in
# the 'MARMICRODB_catalog.filtered.genome.txt'

# path to the big MARMICRODB file:
MARMICRODB_FASTA="MARMICRODB.faa"
# the filtered genome entry IDs:
KEPT_IDS="MARMICRODB_catalog.filtered.genome.txt"

# name of the filtered output:
OUT_FASTA="marmicrodb.filtered.uid.faa"
# this script creates a 'UID2TAX' file for use downstream with the DIAMOND fast protein aligner.
# this is the name of the UID2TAX output:
MARMICRODB_UID2TAX="marmicrodb.filtered.uid2tax.tab"

# run the script:
./process_marmicrodb_fasta.py -k ${KEPT_IDS} -f ${MARMICRODB_FASTA}  -o ${OUT_FASTA} -t ${MARMICRODB_UID2TAX}
# This outputs a filtered fasta file (OUT_FASTA)
# and a uid2tax file for DIAMOND annotation (MARMICRODB_UID2TAX)

# A small number of MARMICRODB.faa sequences contain 
# unexpected numeric values in the amino acid sequence strings:
# example from MARMICRODB.faa:
"	>mmdb1 1A-1_54252
	21826632182675MHGKSGSRGGPVAQPGRALGSHPRGPGFKSRPVHHPLLDALRAWGEL
	>mmdb2 1A-2_54252
	21826632182675MDVLDEVFERVVKARIFRNRSVLSPDYIPDKLPHREREIRALGSIV"

# We employ this script to remove the sequences with numeric values:
MARMICRODB_FASTA="marmicrodb.filtered.uid.faa"
OUTPUT_FASTA="marmicrodb.filtered2.uid.faa"
marmicrodb_remove_numeric_seqs.py -f ${MARMICRODB_FASTA} -o ${OUTPUT_FASTA}

# This removes a total of 323,835 sequences with numeric residues (1.2% of total sequences)
grep -c ">" $MARMICRODB_FASTA # 27890788
grep -c ">" $OUTPUT_FASTA # 27566953

# compress the MMDB uid2tax and FASTA file:
gzip $MARMICRODB_FASTA
gzip $MARMICRODB_UID2TAX
