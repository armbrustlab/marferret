
# download Pfam 34.0 from source:
## Pfam 34.0 (March 2021, 19179 entries) ##
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.full.gz

# Full annotation against Pfam was performed with hmmsearch (HMMER 3.3)
# The protein fasta file for each entry was annotated against Pfam:
for ENTRY in $DATABASE_ENTRIES; do
hmmsearch --cut_tc --domtblout $OUTPUT_DIR/$ENTRY.PFAM34.0.domtblout.tab $PFAM34_PATH $FASTA_PATH
done

# Results are converted to csv format with a custom script:
hmmer_domtblout2csv.py -t $OUTPUT_DIR/$ENTRY.PFAM34.0.domtblout.tab -o $OUTPUT_DIR/$ENTRY.PFAM34.0.domtblout.csv

# Best-scoring Pfam for each frame and contig were retained with a custom script:
"EukRefDB.best_kofam.R"
