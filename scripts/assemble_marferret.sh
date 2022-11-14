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
            printf '\ttranslating nucleotide sequence\n'
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

# Cluster

# filepaths
VERSION=v1
META_FILE=$( realpath "../data/MarFERReT.${VERSION}.metadata.csv" )
RAW_SEQ_DIR=$( realpath ../data/aa_seqs )
WORK_DIR=../data/taxid_grouped
OUTPUT_DIR=../data/clustered
MIN_SEQ_ID=0.99

# make new working directory if need be
if [ ! -d ${WORK_DIR} ]; then 
    mkdir ${WORK_DIR}
fi
WORK_DIR=$( realpath ${WORK_DIR} )

# make new output directory if need be
if [ ! -d ${OUTPUT_DIR} ]; then 
    mkdir ${OUTPUT_DIR}
fi
OUTPUT_DIR=$( realpath ${OUTPUT_DIR} )

# combine amino acid fasta files by taxID and rename MarFERReT protein IDs
UID2TAXID="MarFERReT.${VERSION}.uid2tax.tab"
UID2DEFLINE="MarFERReT.${VERSION}.uid2def.csv"
# move to working directory
pushd $WORK_DIR
# run python script
../../scripts/uniq_id_and_group_by_taxid.py -t ${UID2TAXID} -c ${UID2DEFLINE} \
    -r ${META_FILE} -d ${RAW_SEQ_DIR}
# the UID2TAXID file for use with diamond
gzip $UID2TAXID 
# move both out to data directory
mv "${UID2DEFLINE}" ../
mv "${UID2TAXID}.gz" ../

# cluster proteins from NCBI tax ids with more than one reference sequence
F_TAX_ID=$(head -n1 ${META_FILE} | tr ',' '\n' | grep -Fxn "tax_id" | cut -f1 -d:)
# iterate through NCBI tax ids to be clustered (more than one reference sequence)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -d ); do
    INPUT_FASTA="${TAXID}.combined.faa"
    # make temporary working directory for taxid
    mkdir -p ${TAXID}/${TAXID}_tmp
    # make combined taxid sequence database
    docker run -w /data -v ${WORK_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
        createdb ${INPUT_FASTA} ${TAXID}/${TAXID}.db
    # cluster sequences from combined taxid sequence database
    docker run -w /data -v ${WORK_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
        linclust ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.db \
        ${TAXID}/${TAXID}_tmp --min-seq-id ${MIN_SEQ_ID}
    # select representative sequence from each sequence cluster
    docker run -w /data -v ${WORK_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
        result2repseq ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.db \
        ${TAXID}/${TAXID}.clusters.rep
    # output representative sequence from each sequence cluster
    docker run -w /data -v ${WORK_DIR}:/data ghcr.io/soedinglab/mmseqs2 \
        result2flat ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.db \
        ${TAXID}/${TAXID}.clusters.rep ${TAXID}/${TAXID}.clustered.faa --use-fasta-header
    # moved clustered sequence result to output directory
    mv ${TAXID}/${TAXID}.clustered.faa ${OUTPUT_DIR}/
    # delete temporary working directory
    rm -rf ${TAXID}
done 
# return to script directory
popd

# combine all clustered protein sequence representatives with unclustered
# protein sequences for complete MarFERReT protein database
MARFERRET_FASTA="../data/MarFERReT.${VERSION}.proteins.faa"
# combine NCBI tax IDs with multiple sequence representatives (clustered)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -d ); do
    cat ${OUTPUT_DIR}/${TAXID}.clustered.faa >> ${MARFERRET_FASTA}
done
# combine NCBI tax IDs with single sequence representatives (unclustered)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -u ); do
    cat ${WORK_DIR}/${TAXID}.combined.faa >> ${MARFERRET_FASTA}
done
