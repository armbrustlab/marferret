#!/bin/bash
# AUTHOR: Stephen Blaskowski

# The purpose of this script is to build a diamond database from the 
# MarFERReT database files.


# exit if the script tries to use undeclared variables
set -o nounset
# exit if any pipe commands fail
set -o pipefail
# exit when a command fails
set -o errexit

# input variables
VERSION=v1
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa.gz"
UID2TAXID="MarFERReT.${VERSION}.taxonomies.tab.gz"
DATA_DIR=$( realpath ../data )
TAX_DIR=${DATA_DIR}/diamond/ncbi

# user selects singularity or docker containerization
CONTAINER=""
while [ "${CONTAINER}" == "" ]; do
    printf "\nPlease select an option:\n\n\t1 - singularity\n\t2 - docker\n\nEnter '1' or '2': "
    read selection 
    if [ "${selection}" == "1" ]; then
        echo "Continuing with Singularity containerized workflow"
        CONTAINER="singularity"
    elif [ "${selection}" == "2" ]; then
        echo "Continuing with Docker containerized workflow"
        CONTAINER="docker"
    else 
        printf "\nInvalid selection\n"
    fi
done

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
pushd ${DATA_DIR}
MARFERRET_DMND=diamond/MarFERReT.v1.dmnd
# assign containerized filepaths for taxon notes and names
TAXONNODES=diamond/ncbi/nodes.dmp
TAXONNAMES=diamond/ncbi/names.dmp

# build database with singularity
if [ "${CONTAINER}" == "singularity" ]; then
    singularity exec "${CONTAINER_DIR}/diamond.sif" makedb \
        --in ${MARFERRET_PROTEINS} --db ${MARFERRET_DMND} \
        --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} \
        --taxonmap ${UID2TAXID}
# build database with docker
elif [ "${CONTAINER}" == "docker" ]; then
    docker run -w /data -v $(pwd):/data buchfink/diamond:version2.0.13 makedb \
        --in ${MARFERRET_PROTEINS} --db ${MARFERRET_DMND} \
        --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} \
        --taxonmap ${UID2TAXID}
else
    echo "Containerization not recognized"
    exit
fi

# finished
popd
echo "MarFERReT diamond database construction complete!"
