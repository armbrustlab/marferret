# Download the FASTA files for each entry to this /source_seqs/ directory
# The /source_seqs/ contents should contain the files listed in the 'source_filename'
# column of the metadata.csv file in the /data/ directory

# FASTA files are only necessary for the entries with a 'Y' in the 'accepted' column;
# they will be ingested when calling the 'assemble_marferret.sh' script. Entries with
# a 'N' in the 'accepted' column are not necessary. 
