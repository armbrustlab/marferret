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
WORKING_DIR=../data
# filepaths (relative to working directory)
PFAM_HMMS="pfam/Pfam-A.hmm"
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa"
DOMTBL="pfam/MarFERReT.${VERSION}.annotations.domtblout.tab"
ANNOTATIONS="MarFERReT.${VERSION}.annotations.csv"

# run hmmsearch
docker run -v "$( realpath ../data )":/data biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1 \
    hmmsearch --cut_tc --domtblout ${DOMTBL} $PFAM_HMMS ${MARFERRET_PROTEINS}

# choose best annotation for each protein
./best_kofam.py ${WORKING_DIR}/${DOMTBL} ${WORKING_DIR}/${ANNOTATIONS}
