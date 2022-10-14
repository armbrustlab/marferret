


#### Assigning unique MarFERReT protein IDs and grouping by NCBI taxID

The sequence identifiers (defline) in source sequence FASTA files come in heterogenous formats, may contain extraneous information, and can be problematic for downstream processes. We standardize the MarFERReT protein IDs by assigning an internal sequence identifier, mft[ref_id]_seq[n], where [ref_id] is the unique MarFERReT entry ID, and [n] is a numeric iterator of the sequences within an entry. Entries are also grouped by their NCBI taxID in advance of the clustering step.


`MARFERRET_DIR="/mnt/nfs/projects/ryan/MarFERReT_v1"`

Move to the directory where grouped and sequence-renamed FASTA files will be deposited:
`cd $RAW_SEQ_DIR/taxid_grouped`

Path to the metadata file with ref_id, aa_fasta, tax_id cols:
`REF_DAT="${MARFERRET_DIR}/data/MarFERReT.entry_paths.v1.csv"`

Path to protein sequences:
`RAW_SEQ_DIR="${MARFERRET_DIR}/data/raw_sequence/"`


Name of output file to be created linking unique sequence ID to taxID
`UID2TAXID="${MARFERRET_DIR}/data/EukRefDB.uid2tax.tab"`

Name of output file to be created linking unique sequence ID to source defline
`UID2DEFLINE="${MARFERRET_DIR}/data/EukRefDB.uid2def.csv"`

Run the script to combine entry FASTA files by taxID, assign each sequence a unique ID in the process:
`uniq_id_and_group_by_taxid.py -t ${UID2TAXID} -c ${UID2DEFLINE} -r ${REF_DAT} -d ${RAW_SEQ_DIR}`

#### Intra-taxID clustering

Reduce redundancy within taxIDs with multiple entries by clustering at the 99% protein identity threshold and retaining cluster representatives.

Operate in the directory for taxid-grouped entries from above;

`cd $MARFERRET_DIR/data/raw_sequence/taxid_grouped/`

A simple text list of the 132 multi-entry taxIDs to cluster:

`TAXID2CLUSTER="$MARFERRET_DIR/data/MarFERReT.multi_entry_taxids.txt"`

Output path for clustered sequences:

`OUTPUT_DIR="$MARFERRET_DIR/data/raw_sequence/taxid_grouped/clustered"`

Define the sequence similarity threshold (0.99 = 99%):

`MIN_SEQ_ID=0.99`

Define a function to run mmseqs and linclust together to cluster the input FASTA and generate an output FASTA with the retained cluster representatives. The primary output here is '${TAXID}.nr.aa.fasta', which is copied to the $OUTPUT_DIR.

`function lincluster_marferret {
INPUT_FASTA="MarFERReT.${TAXID}.combined.aa.fasta"
mkdir ${TAXID}
mkdir ${TAXID}/${TAXID}_tmp
cp $MARFERRET_DIR/data/raw_sequence/taxid_grouped/${INPUT_FASTA} ${TAXID}/
mmseqs createdb ${TAXID}/${INPUT_FASTA} ${TAXID}/${TAXID}.db
mmseqs linclust ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.db ${TAXID}/${TAXID}_tmp --min-seq-id ${MIN_SEQ_ID}
mmseqs result2repseq ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.db ${TAXID}/${TAXID}.clusters.rep
mmseqs result2flat ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.db ${TAXID}/${TAXID}.clusters.rep ${TAXID}/${TAXID}.nr.aa.fasta --use-fasta-header
cp ${TAXID}/${TAXID}.nr.aa.fasta ${OUTPUT_DIR}/
if [[ -e ${OUTPUT_DIR}/${TAXID}.nr.aa.fasta ]]; then
rm -rf /scratch/mmseqs/${TAXID}; fi
}`

Iterate through a list of the multi-entry taxIDs to cluster together:

`for TAXID in $(cat ${TAXID2CLUSTER}); do
echo ${TAXID}
lincluster_marferret
done
`

#### Merge all protein sequences into a single file ####

Concatenate the clustered multi-entry taxIDs and the un-clustered single-entry taxIDs sequences together into one file

The clustered multi-entry taxIDs:
`TAXID2CLUSTER="$MARFERRET_DIR/data/MarFERReT.multi_entry_taxids.txt"`
And the other 385 single-entry taxIDs:
`TAXID_UNCLUST="$MARFERRET_DIR/data/MarFERReT.single_entry_taxids.txt"`

Define the name of the single output file:
MARFERRET_FASTA="${MARFERRET_DIR}/data/MarFERReT.v1.proteins.faa"

`cd $MARFERRET_DIR/data/
for TAXID in $(cat ${TAXID2CLUSTER}); do
cat ${MARFERRET_DIR}/data/raw_sequence/taxid_grouped/clustered/${TAXID}.nr.aa.fasta >> ${MARFERRET_FASTA}; done
for TAXID in $(cat ${TAXID_UNCLUST}); do
cat ${MARFERRET_DIR}/data/raw_sequence/taxid_grouped/MarFERReT.${TAXID}.combined.aa.fasta >> ${MARFERRET_FASTA}; done`

MarFERReT.v1.proteins.faa now contains proteins quences from all MarFERReT entries and can be used for downstream analyses.
