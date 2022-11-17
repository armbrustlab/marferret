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
CORES=32
WORKING_DIR=$( realpath ../data )
# filepaths (relative to working directory)
PFAM_HMMS="${WORKING_DIR}/pfam/Pfam-A.hmm"
MARFERRET_PROTEINS="${WORKING_DIR}/MarFERReT.${VERSION}.proteins.faa"
DOMTBL="${WORKING_DIR}/pfam/MarFERReT.${VERSION}.annotations.domtblout.tab"
ANNOTATIONS="${WORKING_DIR}/MarFERReT.${VERSION}.annotations.csv"
CONTAINER_DIR=$( realpath ${WORKING_DIR}/../container )

# run hmmsearch
singularity exec ${CONTAINER_DIR}/hmmer.sif hmmsearch --cpu ${CORES} --cut_tc \
    --domtblout ${DOMTBL} ${PFAM_HMMS} ${MARFERRET_PROTEINS}

# choose best annotation for each protein
singularity exec ${CONTAINER_DIR}/marferret-py.sif ./best_kofam.py \
    ${DOMTBL} ${ANNOTATIONS}
