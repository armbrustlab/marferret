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


# exit if the script tries to use undeclared variables
set -o nounset
# exit if any pipe commands fail
set -o pipefail
# exit when a command fails
set -o errexit

# input variables
CORES=4
VERSION=v1
MARFERRET_DIR=$( realpath ../ )
DATA_DIR="${MARFERRET_DIR}/data"
# filepaths (relative to working directory)
PFAM_HMMS="pfam/Pfam-A.hmm"
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa"
DOMTBL="pfam/MarFERReT.${VERSION}.pfam.domtblout.tab"
ANNOTATIONS="MarFERReT.${VERSION}.best_pfam_annotations.csv"

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

# gunzip ${MARFERRET_PROTEINS} if need be
if [ ! -e "${DATA_DIR}/${MARFERRET_PROTEINS}" ]; then
    if [ -e "${DATA_DIR}/${MARFERRET_PROTEINS}.gz" ]; then
        gunzip "${DATA_DIR}/${MARFERRET_PROTEINS}.gz"
    else
        echo "Could not find input file ${DATA_DIR}/${MARFERRET_PROTEINS}.gz"
    fi
fi

# gunzip ${PFAM_HMMS} if need be
if [ ! -e "${DATA_DIR}/${PFAM_HMMS}" ]; then
    if [ -e "${DATA_DIR}/${PFAM_HMMS}.gz" ]; then
        gunzip "${DATA_DIR}/${PFAM_HMMS}.gz"
    else
        echo "Could not find input file ${DATA_DIR}/${PFAM_HMMS}.gz"
    fi
fi

# make domtblout file
if [ ! -e "${DATA_DIR}/${DOMTBL}" ]; then
    touch "${DATA_DIR}/${DOMTBL}"
fi

# run pfam annotations with singularity
if [ "${CONTAINER}" == "singularity" ]; then
    # run hmmsearch
    singularity exec --no-home --bind ${MARFERRET_DIR}:/marferret \
        ${CONTAINER_DIR}/hmmer.sif hmmsearch --cpu ${CORES} --cut_tc \
        --domtblout "/marferret/data/${DOMTBL}" \
        "/marferret/data/${PFAM_HMMS}" "/marferret/data/${MARFERRET_PROTEINS}"
    # choose best annotation for each protein
    singularity exec --no-home --bind ${MARFERRET_DIR}:/marferret \
        ${CONTAINER_DIR}/marferret-py.sif \
        "/marferret/scripts/python/best_pfam.py" \
        "/marferret/data/${DOMTBL}" "/marferret/data/${ANNOTATIONS}"
# run pfam annotations with docker
elif [ "${CONTAINER}" == "docker" ]; then
    # run hmmsearch
    docker run -v ${DATA_DIR}:/data biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1 \
        hmmsearch --cpu ${CORES} --cut_tc \
        --domtblout ${DOMTBL} ${PFAM_HMMS} ${MARFERRET_PROTEINS}
    # choose best annotation for each protein
    docker run -w /home -v ${MARFERRET_DIR}:/home marferret-py \
        scripts/python/best_pfam.py data/${DOMTBL} data/${ANNOTATIONS}
else
    echo "Containerization not recognized"
    exit
fi

# gzip input proteins
gzip ${DATA_DIR}/${MARFERRET_PROTEINS}
# gzip outputs
gzip ${DATA_DIR}/${DOMTBL}
gzip ${DATA_DIR}/${ANNOTATIONS}
echo "MarFERReT annotation with pfam complete!"
