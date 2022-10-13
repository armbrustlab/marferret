#### Six-frame translation of nucleotide sequence into protein using transeq

Six frame translation was conducted on all sequences ingested in nucleotide-format. A total of 71 entries were translated from nucleotide to protein; these entries are listed with 'nt' in the 'seq_type' field in the primary metadata file: 
[MarFERReT.entry_metadata.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_metadata.v1.csv)

Compressed (gzipped) files were also unzipped at this step for transeq and downstream steps:
`gunzip *.fasta.gz`

For each nucleotide entry ($ENTRY), six-frame translation using the Standard Genetic Code was performed using transeq. Example loop iterating over entry files:
`for ENTRY in $ENTRY_LIST; do
transeq -auto -sformat pearson -frame 6 -sequence ${ENTRY}.fasta -outseq 6tr/${ENTRY}.6tr.fasta
done`

Example filename is shown above; full filenames used internally are listed in [MarFERReT.entry_paths.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_paths.v1.csv)

