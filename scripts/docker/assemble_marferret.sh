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
#         the five fields 'ref_id', 'marferret_name', 'source_filename', 
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

# input variables
VERSION="v1"
MARFERRET_DIR=$( realpath ../ )
MIN_SEQ_ID=0.99     # sequence identity threshold for amino acid clustering
SOURCE_DIR="${MARFERRET_DIR}/data/source_seqs"
SCRIPTS_DIR="${MARFERRET_DIR}/scripts/core"
META_FILE="${MARFERRET_DIR}/data/MarFERReT.${VERSION}.metadata.csv"

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
F_REF_ID=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "ref_id" | cut -f1 -d:)
F_NAME=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "marferret_name" | cut -f1 -d:)
F_FILE=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "source_filename" | cut -f1 -d:)
F_FASTA=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "aa_fasta" | cut -f1 -d:)
F_SEQ_TYPE=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "seq_type" | cut -f1 -d:)

# iterate through nucleotide sequences listed in metatdata file
while IFS=',' read -r ref_id marferret_name source_filename seq_type aa_fasta; do 
    # check that the sequence is in the source_seqs directory
    if [ ! -e "${SOURCE_DIR}/${source_filename}" ]; then 
        echo "WARNING: Filename ${source_filename} not found."
        echo "See missing_source_seqs.txt for complete list of missing sequences."
        echo "${source_filename} not found in ${SOURCE_DIR}" >> missing_source_seqs.txt
    else
        echo $ref_id $marferret_name $source_filename $seq_type
        # build standardized sequence name
        seq_name="${aa_fasta%.*}"
        # rename amino acids and move to amino acid sequence directory
        if [ "${seq_type}" == "aa" ]; then
            cat "${SOURCE_DIR}/${source_filename}" >> "${AA_DIR}/${aa_fasta}"
        elif [ "${seq_type}" == "nt" ]; then
            printf '\ttranslating nucleotide sequence\n'
            # run transeq to translate nucleotide sequences into proteins
            docker run -v ${MARFERRET_DIR}:/data \
                biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1 \
                transeq -auto -sformat pearson -frame 6 \
                -sequence "data/source_seqs/${source_filename}" \
                -outseq "data/temp/${seq_name}.6tr.faa"
            # run frame selection to select longest coding frame
            docker run -w /home -v ${MARFERRET_DIR}:/home marferret-py \
                ${SCRIPTS_DIR}/keep_longest_frame.py -l 1 \
                "data/temp/${seq_name}.6tr.faa"
            # move frame selected file to amino acid sequence directory
            mv "${TMP_DIR}/${seq_name}.6tr.bf1.faa" "${AA_DIR}/${aa_fasta}"
        fi
    fi
done < <( tail -n +2 ${META_FILE} | cut -f"${F_REF_ID},${F_NAME},${F_SEQ_TYPE},${F_FILE},${F_FASTA}" -d, )

# clean up temp directory
rm -rf ${TMP_DIR}

# combine amino acid fasta files by taxID and rename MarFERReT protein IDs
# make new directory for taxid combined sequences
TAX_DIR="${MARFERRET_DIR}/data/taxid_grouped"
if [ ! -d ${TAX_DIR} ]; then 
    mkdir ${TAX_DIR}
fi
# run python script
docker run -w /home -v ${MARFERRET_DIR}:/home marferret-py \
    ${SCRIPTS_DIR}/group_by_taxid.py data/aa_seqs \
    "data/MarFERReT.${VERSION}.metadata.csv" \
    -o data/taxid_grouped
# move mapping files to the data directory
mv "${TAX_DIR}/taxonomies.tab" "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.taxonomies.tab"
mv "${TAX_DIR}/proteins_info.tab" "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.proteins_info.tab"
# gzip the taxonomies.tab file for use with diamond
gzip "${MARFERRET_DIR}/data/MarFERReT.${VERSION}.taxonomies.tab"

# cluster proteins from NCBI tax ids with more than one reference sequence
F_TAX_ID=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "tax_id" | cut -f1 -d:)
# make new directory for clustered sequences
CLUSTER_DIR="${MARFERRET_DIR}/data/clustered"
if [ ! -d ${CLUSTER_DIR} ]; then 
    mkdir ${CLUSTER_DIR}
fi
# move to work in TAX_DIR directory
pushd ${TAX_DIR}
# iterate through NCBI tax ids to be clustered (more than one reference sequence)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -d ); do
    INPUT_FASTA="${TAXID}.combined.faa"
    # make temporary working directory for taxid
    mkdir -p ${TAXID}/${TAXID}_tmp
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
        ${TAXID}/${TAXID}.clusters.rep ${TAXID}/${TAXID}.clustered.faa --use-fasta-header
    # moved clustered sequence result to output directory
    mv ${TAXID}/${TAXID}.clustered.faa ${CLUSTER_DIR}/
    # delete temporary working directory
    rm -rf ${TAXID}
done 
# return to original directory
popd

# combine all clustered protein sequence representatives with unclustered
# protein sequences for complete MarFERReT protein database
MARFERRET_FASTA="../data/MarFERReT.${VERSION}.proteins.faa"
# combine NCBI tax IDs with multiple sequence representatives (clustered)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -d ); do
    cat ${CLUSTER_DIR}/${TAXID}.clustered.faa >> ${MARFERRET_FASTA}
done
# combine NCBI tax IDs with single sequence representatives (unclustered)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -u ); do
    cat ${TAX_DIR}/${TAXID}.combined.faa >> ${MARFERRET_FASTA}
done
