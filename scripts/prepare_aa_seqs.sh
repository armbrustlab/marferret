#!/bin/bash
# AUTHOR: Stephen Blaskowski

# This script performs six-frame translation of nucleotide sequence into 
# amino acid sequences using transeq, and then selects the longest uninterrupted
# amino acid sequence for futher analysis using a custom python script.

# Six frame translation was conducted on all sequences ingested in nucleotide 
# format. A total of 71 entries were translated from nucleotide to protein; 
# these entries are listed with 'nt' in the 'seq_type' field in the primary 
# metadata file: MarFERReT.entry_metadata.v1.csv

# Place all nucleotide sequences in the nucleotide directory listed here. Final 
# frame selected amino acid sequences will be output to the amino acid sequence 
# directory listed here. Note that any compressed files (e.g. gzipped) should
# be uncompressed prior to running this script.

VERSION=v1
SOURCE_DIR=$( realpath ../data/source_seqs )
META_FILE=$( realpath "../data/MarFERReT.${VERSION}.metadata.csv" )

# make temp directory to put 6 frame translated files
TMP_DIR=../data/temp
if [ ! -d "${TMP_DIR}" ]; then
    mkdir "${TMP_DIR}"
fi
TMP_DIR=$( realpath ${TMP_DIR} )

# make directory for amino acids
AA_DIR=../data/aa_seqs
if [ ! -d "${AA_DIR}" ]; then
    mkdir "${AA_DIR}"
fi
AA_DIR=$( realpath ${AA_DIR} )

# get field numbers from metadata file
F_REF_ID=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "ref_id" | cut -f1 -d:)
F_NAME=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "marferret_name" | cut -f1 -d:)
F_FILE=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "source_filename" | cut -f1 -d:)
F_SEQ_TYPE=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "seq_type" | cut -f1 -d:)

# iterate through nucleotide sequences listed in metatdata file
while IFS=',' read -r ref_id marferret_name source_filename seq_type; do 
    # check that the sequence is in the source_seqs directory
    if [ ! -e "${SOURCE_DIR}/${source_filename}" ]; then 
        echo "WARNING: Filename ${source_filename} not found."
        echo "See missing_source_seqs.txt for complete list of missing sequences."
        echo "${source_filename} not found in ${SOURCE_DIR}" >> missing_source_seqs.txt
    else
        echo $ref_id $marferret_name $source_filename $seq_type
        # build standardized sequence name
        seq_name="${ref_id}_${marferret_name}"
        # rename amino acids and move to amino acid sequence directory
        if [ "${seq_type}" == "aa" ]; then
            cat "${SOURCE_DIR}/${source_filename}" >> "${AA_DIR}/${seq_name}.faa"
        elif [ "${seq_type}" == "nt" ]; then
            printf '\ttranslating nucleotide sequence'
            # run transeq to translate nucleotide sequences into proteins
            docker run -v "$( realpath ../data )":/data \
                biocontainers/emboss:v6.6.0dfsg-7b1-deb_cv1 \
                transeq -auto -sformat pearson -frame 6 \
                -sequence "/data/source_seqs/${source_filename}" \
                -outseq "/data/temp/${seq_name}.6tr.faa"
            # run frame selection to select longest coding frame
            ./keep_longest_frame.py -l 1 "${TMP_DIR}/${seq_name}.6tr.faa"
            # move frame selected file to amino acid sequence directory
            mv "${TMP_DIR}/${seq_name}.6tr.bf1.faa" "${AA_DIR}/${seq_name}.faa"
        fi
    fi
done < <( tail -n +2 ${META_FILE} | cut -f"${F_REF_ID},${F_NAME},${F_SEQ_TYPE},${F_FILE}" -d, )

# clean up temp directory
rm -rf ${TMP_DIR}
