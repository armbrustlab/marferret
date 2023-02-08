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


# # exit if the script tries to use undeclared variables
# set -o nounset
# # exit if any pipe commands fail
# set -o pipefail
# # exit when a command fails
# set -o errexit

# input variables
VERSION=v1
WORKING_DIR=$( realpath ../data )
# filepaths (relative to working directory)
PFAM_HMMS="pfam/Pfam-A.hmm"
MARFERRET_PROTEINS="MarFERReT.${VERSION}.proteins.faa.gz"
DOMTBL="pfam/MarFERReT.${VERSION}.pfam.domtblout.tab"
ANNOTATIONS="MarFERReT.${VERSION}.pfam.csv"

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

# run pfam annotations with singularity
if [ "${CONTAINER}" == "singularity" ]; then
    # run hmmsearch
    singularity exec ${CONTAINER_DIR}/hmmer.sif hmmsearch --cpu ${CORES} --cut_tc \
        --domtblout ${DOMTBL} ${PFAM_HMMS} ${MARFERRET_PROTEINS}
    # choose best annotation for each protein
    singularity exec ${CONTAINER_DIR}/marferret-py.sif ./best_kofam.py \
        ${DOMTBL} ${ANNOTATIONS}

# run pfam annotations with docker
elif [ "${CONTAINER}" == "docker" ]; then
    # run hmmsearch
    docker run -v ${WORKING_DIR}:/data biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1 \
        hmmsearch --cut_tc --domtblout ${DOMTBL} ${PFAM_HMMS} ${MARFERRET_PROTEINS}
    # choose best annotation for each protein
    docker run -w /home -v ${WORKING_DIR}:/home marferret-py ../scripts/python/best_kofam.py \
        ${DOMTBL} ${ANNOTATIONS}
else
    echo "Containerization not recognized"
    exit
fi

# gzip outputs
gzip ${DOMTBL}
gzip ${ANNOTATIONS}
echo "MarFERReT annotation with pfam complete!"
