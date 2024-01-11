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

# directories
MARFERRET_DIR=$( realpath ../ )
DATA_DIR="${MARFERRET_DIR}/data"
TAX_DIR="${DATA_DIR}/diamond/ncbi"

# pull version from metadata file
FILENAME=$( basename ${DATA_DIR}/MarFERReT.*.metadata.csv )
FRONT="${FILENAME%%.metadata.csv}"
VERSION="${FRONT##MarFERReT.}"

# input files
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa.gz"
UID2TAXID="MarFERReT.${VERSION}.taxonomies.tab.gz"

# user selects singularity or docker containerization
CONTAINER=""
while [ "${CONTAINER}" == "" ]; do
    printf "\nPlease select an option:\n\n\t1 - singularity\n\t2 - docker\n\nEnter '1' or '2': "
    read selection 
    if [ "${selection}" == "1" ]; then
        echo "Continuing with Singularity containerized workflow"
        CONTAINER="singularity"
        CONTAINER_DIR="${MARFERRET_DIR}/containers"
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
if [ ! -e "taxdmp.zip" ]; then
    wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
fi
if [ ! -e "nodes.dmp" ]; then
    unzip taxdmp.zip
fi
popd

# build diamond database
pushd ${DATA_DIR}
MARFERRET_DMND=diamond/MarFERReT.v1.dmnd
# assign containerized filepaths for taxon notes and names
TAXONNODES=diamond/ncbi/nodes.dmp
TAXONNAMES=diamond/ncbi/names.dmp

# build database with singularity
if [ "${CONTAINER}" == "singularity" ]; then
    singularity exec --no-home --bind ${DATA_DIR} \
        "${CONTAINER_DIR}/diamond.sif" diamond makedb \
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
