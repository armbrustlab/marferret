#!/bin/bash
# AUTHOR: Stephen Blaskowski

# The purpose of this script is to build a diamond database from the 
# MarFERReT database files.


VERSION=v1
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa"
UID2TAXID="MarFERReT.${VERSION}.uid2tax.tab.gz"
TAX_DIR=../data/diamond/ncbi

# make directory for diamond database and ncbi tax data
if [ ! -d "${TAX_DIR}" ]; then
    mkdir -p "${TAX_DIR}"
fi

# download and unzip taxdmp.zip files from NCBI
pushd ${TAX_DIR}
wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
unzip taxdmp.zip
popd

# build diamond database
pushd ../data
MARFERRET_DMND=diamond/MarFERReT.v1.dmnd
# assign containerized filepaths for taxon notes and names
TAXONNODES=diamond/ncbi/nodes.dmp
TAXONNAMES=diamond/ncbi/names.dmp
# build database
docker run -w /data -v $(pwd):/data buchfink/diamond:version2.0.13 makedb \
    --in ${MARFERRET_PROTEINS} --db ${MARFERRET_DMND} \
    --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} \
    --taxonmap ${UID2TAXID}
popd
