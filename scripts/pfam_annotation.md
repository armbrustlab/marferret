#### Pfam annotation of the reference database

Work in the pfam directory for MarFERReT:
`cd ${MARFERRET_DIR}/pfam`

MarFERReT v1 was built using Pfam 34.0 (downloaded March 2021, 19,179 total hmm profiles)
Download Pfam hmm profiles from source:
`wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz
gunzip Pfam-A.full.gz`


Full annotation against Pfam is performed with hmmsearch (HMMER 3.3)
Path to Pfam hmm profiles. For bulk annotation using hmmsearch, the hmm file works best unzipped:
`PFAM_HMMS="${MARFERRET_DIR}/pfam/Pfam-A.hmm"`

The protein fasta file for each entry was annotated against Pfam using this code, using the internal trusted-cutoff score (--cut_tc) as a minimum threshold:
`# $MARFERRET_PROTEINS can be a simple text file with a list of entry handles and paths
# Example loop:
for ENTRY in $MARFERRET_PROTEINS; do
hmmsearch --cut_tc --domtblout $OUTPUT_DIR/$ENTRY.PFAM34.0.domtblout.tab $PFAM_HMMS ${ENTRY}.aa.fasta
done`

Results are converted to csv format with a custom script:
hmmer_domtblout2csv.py -t $OUTPUT_DIR/$ENTRY.PFAM34.0.domtblout.tab -o $OUTPUT_DIR/$ENTRY.PFAM34.0.domtblout.csv

Best-scoring Pfam for each frame and contig were retained with a custom script:
"EukRefDB.best_kofam.R"
