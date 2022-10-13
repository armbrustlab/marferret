#### Six-frame translation of nucleotide sequence into protein using transeq

Six frame translation was conducted on all sequences ingested in nucleotide-format. A total of 71 entries were translated from nucleotide to protein; these entries are listed with 'nt' in the 'seq_type' field in the primary metadata file: 
[MarFERReT.entry_metadata.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_metadata.v1.csv)


For each nucleotide entry ($ENTRY), six-frame translation using the Standard Genetic Code was performed using transeq: 
`transeq -auto -sformat pearson -frame 6 -sequence ${ENTRY}.fasta -outseq 6tr/${ENTRY}.6tr.fasta`

Exampe filename is shown above; full filenames used internally are listed in [MarFERReT.entry_paths.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_paths.v1.csv)

