#!/bin/bash
# AUTHOR: Stephen Blaskowski

# The purpose of this script is to cluster amino acid sequences for reference
# genomes with the same NCBI taxID. 

# During the course of this, the script will also rename protein sequences to 
# standardize protein IDs. We standardize the MarFERReT protein IDs by 
# assigning an internal sequence identifier, mft[ref_id]_seq[n], where [ref_id] 
# is the unique MarFERReT entry ID, and [n] is a numeric iterator of the 
# sequences within an entry. Entries are then grouped by their NCBI taxID in 
# advance of the clustering step.


# filepaths
META_FILE=$( realpath ../data/MarFERReT.metadata.v1.csv )
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
UID2TAXID='MarFERReT.v1.uid2tax.tab'
UID2DEFLINE='MarFERReT.v1.uid2def.csv'
# move to working directory
pushd $WORK_DIR
# run python script
../scripts/uniq_id_and_group_by_taxid.py -t ${UID2TAXID} -c ${UID2DEFLINE} \
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
    INPUT_FASTA="MarFERReT.${TAXID}.combined.aa.fasta"
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
        ${TAXID}/${TAXID}.clusters.rep ${TAXID}/${TAXID}.nr.aa.fasta --use-fasta-header
    # moved clustered sequence result to output directory
    mv ${TAXID}/${TAXID}.nr.aa.fasta ${OUTPUT_DIR}/
    # delete temporary working directory
    rm -rf ${TAXID}
done 
# return to script directory
popd

# combine all clustered protein sequence representatives with unclustered
# protein sequences for complete MarFERReT protein database
MARFERRET_FASTA="../data/MarFERReT.v1.proteins.faa"
# combine NCBI tax IDs with multiple sequence representatives (clustered)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -d ); do
    cat ${OUTPUT_DIR}/${TAXID}.nr.aa.fasta >> ${MARFERRET_FASTA}
done
# combine NCBI tax IDs with single sequence representatives (unclustered)
for TAXID in $( tail -n +2 $META_FILE | cut -d, -f $F_TAX_ID | sort | uniq -u ); do
    cat ${WORK_DIR}/MarFERReT.${TAXID}.combined.aa.fasta >> ${MARFERRET_FASTA}
done

