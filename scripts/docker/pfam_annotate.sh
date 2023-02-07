#!/bin/bash
# AUTHOR: Stephen Blaskowski

# The purpose of this script is to annotate the database proteins

# If you do not have the latest Pfam database downloaded, please make a pfam
# directory under the data directory (../data/pfam) and download the pfam 
# database here by going to the website `http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/` 
# or by using the following commands:

# wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 

# Use the md5 checksum values at the ftp site to ensure that the full
# database was downloaded successfully.

# wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/md5_checksums
# md5sum --ignore-missing -c md5_checksums

# Make sure the hammer profiles are unzipped prior to using hmmsearch
# gunzip Pfam-A.hmm.gz


VERSION=v1
WORKING_DIR=$( realpath ../../data )
# filepaths (relative to working directory)
PFAM_HMMS="pfam/Pfam-A.hmm"
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa"
DOMTBL="pfam/MarFERReT.${VERSION}.pfam.domtblout.tab"
ANNOTATIONS="MarFERReT.${VERSION}.pfam.csv"

# run hmmsearch
docker run -v ${WORKING_DIR}:/data biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1 \
    hmmsearch --cut_tc --domtblout ${DOMTBL} ${PFAM_HMMS} ${MARFERRET_PROTEINS}

# choose best annotation for each protein
docker run -w /home -v ${WORKING_DIR}:/home marferret-py ../core/best_kofam.py \
    ${DOMTBL} ${ANNOTATIONS}
