#!/bin/bash
# AUTHOR: Stephen Blaskowski

# This script assembles the MarFERReT database of reference protein sequences
# from input reference genomes.

# Dependencies: 
#   Ensure that the software containers needed to run the assembly steps have 
#   been set up on your machine (see containers/build_docker_images.sh)

# Inputs: 
#   - Reference genomes in the data/source_seqs/ directory
#   - data/MarFERReT.entry_metadata.v1.csv
#       * NOTE: The MarFERReT.entry_metadata.v1.csv file must contain at least
#         the five fields 'entry_id', 'marferret_name', 'source_filename', 
#         'aa_fasta' and 'seq_type', and the filename listed under 
#         'source_filename' must match the filename of the reference genome 
#         in data/source_seqs/ 

# Outputs:
#   - data/MarFERReT.${VERSION}.proteins.faa
#       * Clustered MarFERReT protein database
#   - data/MarFERReT.${VERSION}.taxonomies.tab.gz
#   - data/MarFERReT.${VERSION}.proteins_info.tab
#   - data/aa_seqs directory with translated & standardized amino acid sequences
#   - data/taxid_grouped directory with amino acid sequences grouped by taxid
#   - data/clustered directory with amino acid sequences clustered within taxid

# This script first performs six-frame translation of all input nucleotide 
# sequence into amino acid sequences using transeq, and then selects the 
# longest uninterrupted amino acid sequence for futher analysis using a custom 
# python script. These sequences are standardized and combined with the input
# amino acid sequences in the `aa_seqs` directory. Next, all amino acid 
# sequences are renamed with a custom MarFERReT unique identifier, which can 
# be mapped to its original name, as well as the reference sequence's 
# taxonomical identifier in the `taxonomies.tab` and `proteins_info.tab` files. 
# All amino acid sequences corresponding to the same taxonomical id are then 
# grouped together into a fasta file in the `taxid_grouped` directory. The 
# script then uses mmseqs2 to cluster amino acids within a taxid group and 
# output the non-redundant cluster representatives to fasta files in the 
# `clustered` directory. Finally, all MarFERReT protein sequences are
# concatenated into the MarFERReT.${VERSION}.proteins.faa file.

# exit if the script tries to use undeclared variables
set -o nounset
# exit if any pipe commands fail
set -o pipefail
# exit when a command fails
set -o errexit

# input variables
VERSION="v1"
MARFERRET_DIR=$( realpath ../ )
MIN_SEQ_ID=0.99     # sequence identity threshold for amino acid clustering
SOURCE_DIR="${MARFERRET_DIR}/data/source_seqs"
META_FILE="${MARFERRET_DIR}/data/MarFERReT.${VERSION}.metadata.csv"

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

# make directory for amino acids
AA_DIR="${MARFERRET_DIR}/data/aa_seqs"
if [ ! -d "${AA_DIR}" ]; then
    mkdir "${AA_DIR}"
fi

# make temp working directory
TMP_DIR="${MARFERRET_DIR}/data/temp"
if [ ! -d "${TMP_DIR}" ]; then
    mkdir "${TMP_DIR}"
fi

# get field numbers from metadata file
F_ENTRY_ID=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "entry_id" | cut -f1 -d:)
F_NAME=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "marferret_name" | cut -f1 -d:)
F_TAX_ID=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "tax_id" | cut -f1 -d:)
F_FILE=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "source_filename" | cut -f1 -d:)
F_SEQ_TYPE=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "seq_type" | cut -f1 -d:)
F_FASTA=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "aa_fasta" | cut -f1 -d:)

# iterate through nucleotide sequences listed in metatdata file
while IFS=',' read -r entry_id marferret_name source_filename seq_type aa_fasta; do 
    # check that the sequence is in the source_seqs directory
    if [ ! -e "${SOURCE_DIR}/${source_filename}" ]; then 
        echo "WARNING: Filename ${source_filename} not found."
        echo "See missing_source_seqs.txt for complete list of missing sequences."
        echo "${source_filename} not found in ${SOURCE_DIR}" >> missing_source_seqs.txt
    else
        echo $entry_id $marferret_name $source_filename $seq_type
        # build standardized sequence name
        seq_name="${aa_fasta%.*}"
        # rename amino acids and move to amino acid sequence directory
        if [ "${seq_type}" == "aa" ]; then
            cat "${SOURCE_DIR}/${source_filename}" >> "${AA_DIR}/${aa_fasta}"
        elif [ "${seq_type}" == "nt" ]; then
            printf '\ttranslating nucleotide sequence\n'
            # run transeq to translate nucleotide sequences into proteins
            if [ "${CONTAINER}" == "singularity" ]; then
                singularity exec --no-home --bind ${MARFERRET_DIR}:/marferret \
                    "${CONTAINER_DIR}/emboss.sif" \
                    transeq -auto -sformat pearson -frame 6 \
                    -sequence "/marferret/data/source_seqs/${source_filename}" \
                    -outseq "/marferret/data/temp/${seq_name}.6tr.faa"
            elif [ "${CONTAINER}" == "docker" ]; then
                docker run -v ${MARFERRET_DIR}:/data \
                    biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1 \
                    transeq -auto -sformat pearson -frame 6 \
                    -sequence "data/source_seqs/${source_filename}" \
                    -outseq "data/temp/${seq_name}.6tr.faa"
            else
                echo "Containerization not recognized"
                exit
            fi
            # run frame selection to select longest coding frame
            if [ "${CONTAINER}" == "singularity" ]; then
                singularity exec --no-home --bind ${MARFERRET_DIR}:/marferret \
                    "${CONTAINER_DIR}/marferret-py.sif" \
                    "/marferret/scripts/python/keep_longest_frame.py" -l 1 \
                    "/marferret/data/temp/${seq_name}.6tr.faa"
            elif [ "${CONTAINER}" == "docker" ]; then
                docker run -w /home -v ${MARFERRET_DIR}:/home marferret-py \
                    scripts/python/keep_longest_frame.py -l 1 \
                    "data/temp/${seq_name}.6tr.faa"
            else
                echo "Containerization not recognized"
                exit
            fi
            # move frame selected file to amino acid sequence directory
            mv "${TMP_DIR}/${seq_name}.6tr.bf1.faa" "${AA_DIR}/${aa_fasta}"
        fi
    fi
done < <( tail -n +2 ${META_FILE} | cut -f"${F_ENTRY_ID},${F_NAME},${F_SEQ_TYPE},${F_FILE},${F_FASTA}" -d, )
# clean up temp directory
rm -rf ${TMP_DIR}

# combine amino acid fasta files by taxID and rename MarFERReT protein IDs
# make new directory for taxid combined sequences
TAX_DIR="${MARFERRET_DIR}/data/taxid_grouped"
if [ ! -d ${TAX_DIR} ]; then 
    mkdir ${TAX_DIR}
fi
# run python script group_by_taxid.py
if [ "${CONTAINER}" == "singularity" ]; then
    singularity exec --no-home --bind ${MARFERRET_DIR}:/marferret \
        "${CONTAINER_DIR}/marferret-py.sif" \
        "/marferret/scripts/python/group_by_taxid.py" \
        "/marferret/data/aa_seqs" \
        "/marferret/data/MarFERReT.${VERSION}.metadata.csv" \
        -o "/marferret/data/taxid_grouped"
elif [ "${CONTAINER}" == "docker" ]; then
    docker run -w /home -v ${MARFERRET_DIR}:/home marferret-py \
        scripts/python/group_by_taxid.py data/aa_seqs \
        "data/MarFERReT.${VERSION}.metadata.csv" \
        -o data/taxid_grouped
else
    echo "Containerization not recognized"
    exit
fi
# move mapping files to the data directory
mv "${TAX_DIR}/taxonomies.tab" "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.taxonomies.tab"
mv "${TAX_DIR}/proteins_info.tab" "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.proteins_info.tab"
# gzip mapping files
gzip "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.taxonomies.tab"
gzip "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.proteins_info.tab"

# cluster proteins from NCBI tax ids with more than one reference sequence
# make new directory for clustered sequences
CLUSTER_DIR="${MARFERRET_DIR}/data/clustered"
if [ ! -d ${CLUSTER_DIR} ]; then 
    mkdir ${CLUSTER_DIR}
fi
# move to work in TAX_DIR directory
pushd ${TAX_DIR}
# iterate through unique NCBI tax ids
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq ); do
    INPUT_FASTA="${TAXID}.combined.faa"
    # make temporary working directory for taxid
    mkdir -p ${TAXID}/${TAXID}_tmp
    # run clustering in singularity
    if [ "${CONTAINER}" == "singularity" ]; then
        # make combined taxid sequence database
        singularity exec --no-home --bind ${TAX_DIR}:/tax_dir \
            "${CONTAINER_DIR}/mmseqs2.sif" mmseqs_avx2 createdb \
            "/tax_dir/${INPUT_FASTA}" "/tax_dir/${TAXID}/${TAXID}.db"
        # cluster sequences from combined taxid sequence database
        singularity exec --no-home --bind ${TAX_DIR}:/tax_dir \
            "${CONTAINER_DIR}/mmseqs2.sif" mmseqs_avx2 linclust \
            "/tax_dir/${TAXID}/${TAXID}.db" \
            "/tax_dir/${TAXID}/${TAXID}.clusters.db" \
            "/tax_dir/${TAXID}/${TAXID}_tmp" --min-seq-id ${MIN_SEQ_ID}
        # select representative sequence from each sequence cluster
        singularity exec --no-home --bind ${TAX_DIR}:/tax_dir \
            "${CONTAINER_DIR}/mmseqs2.sif" mmseqs_avx2 result2repseq \
            "/tax_dir/${TAXID}/${TAXID}.db" \
            "/tax_dir/${TAXID}/${TAXID}.clusters.db" \
            "/tax_dir/${TAXID}/${TAXID}.clusters.rep"
        # output representative sequence from each sequence cluster
        singularity exec --no-home --bind ${TAX_DIR}:/tax_dir \
            "${CONTAINER_DIR}/mmseqs2.sif" mmseqs_avx2 result2flat \
            "/tax_dir/${TAXID}/${TAXID}.db" "/tax_dir/${TAXID}/${TAXID}.db" \
            "/tax_dir/${TAXID}/${TAXID}.clusters.rep" \
            "/tax_dir/${TAXID}/${TAXID}.clustered.faa" --use-fasta-header
    # run clustering in docker
    elif [ "${CONTAINER}" == "docker" ]; then
        # make combined taxid sequence database
        docker run -w /data -v ${TAX_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
            createdb ${INPUT_FASTA} ${TAXID}/${TAXID}.db
        # cluster sequences from combined taxid sequence database
        docker run -w /data -v ${TAX_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
            linclust ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.db \
            ${TAXID}/${TAXID}_tmp --min-seq-id ${MIN_SEQ_ID}
        # select representative sequence from each sequence cluster
        docker run -w /data -v ${TAX_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
            result2repseq ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.db \
            ${TAXID}/${TAXID}.clusters.rep
        # output representative sequence from each sequence cluster
        docker run -w /data -v ${TAX_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
            result2flat ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.db \
            ${TAXID}/${TAXID}.clusters.rep \
            ${TAXID}/${TAXID}.clustered.faa --use-fasta-header
    else
        echo "Containerization not recognized"
        exit
    fi
    # moved clustered sequence result to output directory
    mv ${TAXID}/${TAXID}.clustered.faa ${CLUSTER_DIR}/
    # delete temporary working directory
    rm -rf ${TAXID}
done 
# return to original directory
popd

# combine all clustered protein sequence representatives to make completed
# MarFERReT protein database
MARFERRET_FASTA="${MARFERRET_DIR}/data/MarFERReT.${VERSION}.proteins.faa"
# combine clustered NCBI tax IDs 
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq ); do
    cat ${CLUSTER_DIR}/${TAXID}.clustered.faa >> ${MARFERRET_FASTA}
done

# gzip output MarFERReT.${VERSION}.proteins.faa file
gzip ${MARFERRET_FASTA}
echo "MarFERReT database construction complete!"
