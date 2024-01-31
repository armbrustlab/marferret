#!/bin/bash

# AUTHOR: Ryan Groussman, PhD

# Navigate to the marmicrodb directory in the local MARFERRET_DIR
cd ${MARFERRET_DIR}/data/marmicrodb

# Define the location of the MarFERReT and MARMICRODB protein FASTA files:
MARFERRET_SEQS="${MARFERRET_DIR}/data/MarFERReT.v1.proteins.faa.gz"
MARMICRODB_SEQS="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered2.uid.faa.gz"

# Approximate sizes of each:
ls -lht $MARFERRET_SEQS # 4.4G
ls -lht $MARMICRODB_SEQS # 5.0G

# Simple concatenation of the gzipped files:
cat $MARFERRET_SEQS $MARMICRODB_SEQS > MarFERReT.MARMICRODB.v1.1.combined.faa.gz
COMBINED_SEQS="${MARFERRET_DIR}/data/marmicrodb/MarFERReT.MARMICRODB.v1.1.combined.faa.gz"
ls -lht $COMBINED_SEQS # 9.4G
zgrep -c ">" $COMBINED_SEQS # 55517966

# Define paths of the taxonomy table files:
# the order matters here; the MarFERReT.v1.taxonomies.tab.gz contains
# the necessary header, and the marmicrodb.filtered.uid2tax.tab.gz
# file does not contain a header so that it can be appended.
MARFERRET_UID2TAX="${MARFERRET_DIR}/data/MarFERReT.v1.taxonomies.tab.gz"
MARMICRODB_UID2TAX="${MARFERRET_DIR}/data/marmicrodb/marmicrodb.filtered.uid2tax.tab.gz"

# Simple concatenation of the gzipped files:
cat $MARFERRET_UID2TAX $MARMICRODB_UID2TAX > MarFERReT.MARMICRODB.v1.1.combined.uid2tax.tab.gz

# approximate size:
ls -lht MarFERReT.MARMICRODB.v1.1.combined.uid2tax.tab.gz
# 156M


#### Build the DIAMOND database
# Navigate the the DIAMOND directory:
cd ${MARFERRET_DIR}/data/marmicrodb/dmnd/
NCORES=4 # adjust this value for your specific system


# Define locations of combined FASTA and table files:
COMBINED_SEQS="${MARFERRET_DIR}/data/marmicrodb/MarFERReT.MARMICRODB.v1.1.combined.faa.gz"
UID2TAXID="${MARFERRET_DIR}/data/marmicrodb/MarFERReT.MARMICRODB.v1.1.combined.uid2tax.tab.gz"

# Locations of NCBI Taxonomy files:
TAXONNODES="${MARFERRET_DIR}/data/diamond/ncbi/nodes.dmp"
TAXONNAMES="${MARFERRET_DIR}/data/diamond/ncbi/names.dmp"
MFT_MMDB_DMND_DB="${MARFERRET_DIR}/data/marmicrodb/dmnd/MarFERReT.MARMICRODB.v1.1.combined.dmnd"

DATA_DIR=${MARFERRET_DIR}/data
TAX_DIR=${DATA_DIR}/diamond/ncbi
CONTAINER_DIR="${MARFERRET_DIR}/containers"


singularity exec --no-home --bind ${DATA_DIR} \
        "${CONTAINER_DIR}/diamond.sif" diamond makedb \
        --in ${COMBINED_SEQS} --db ${MFT_MMDB_DMND_DB} \
        --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} \
        --taxonmap ${UID2TAXID}

