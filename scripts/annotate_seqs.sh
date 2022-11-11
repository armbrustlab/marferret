#!/bin/bash
# AUTHOR: Stephen Blaskowski

# The purpose of this script is to annotate the database proteins

# If you do not have the latest Pfam database downloaded, please make a pfam
# directory under the data directory (../data/pfam) and download the pfam 
# database here by going to the website `http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/` 
# or by using the following commands:

# wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
# gunzip Pfam-A.hmm.gz

# Use the md5 checksum values at the ftp site to ensure that the full
# database was downloaded successfully.


VERSION=v1
PFAM_HMMS=../data/pfam/Pfam-A.hmm
MARFERRET_PROTEINS="../data/MarFERReT.${VERSION}.proteins.faa"
OUTPUT_DIR=../data/pfam

# run hmmsearch
DOMTBL="$OUTPUT_DIR/MarFERReT.${VERSION}.annotations.domtblout.tab"
hmmsearch --cut_tc --domtblout ${DOMTBL} $PFAM_HMMS ${MARFERRET_PROTEINS}

# choose best annotation for each protein
ANNOTATIONS="../data/MarFERReT.${VERSION}.annotations.csv"
./best_kofam.py ${DOMTBL} ${ANNOTATIONS}
