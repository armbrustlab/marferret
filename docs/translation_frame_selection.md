#### Six-frame translation of nucleotide sequence into protein using transeq

Six frame translation was conducted on all sequences ingested in nucleotide-format. A total of 71 entries were translated from nucleotide to protein; these entries are listed with 'nt' in the 'seq_type' field in the primary metadata file: 
[MarFERReT.entry_metadata.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_metadata.v1.csv)

Links to the nucleotide files: [MarFERReT.entry_paths.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_paths.v1.csv)

Compressed (gzipped) files were also unzipped at this step for transeq and downstream steps:
`gunzip *.fasta.gz`

For each nucleotide entry ($ENTRY), six-frame translation using the Standard Genetic Code was performed using transeq. Example loop iterating over entry files:
```
for ENTRY in $ENTRY_LIST; do
# run transeq:
transeq -auto -sformat pearson -frame 6 -sequence ${ENTRY}.fasta -outseq 6tr/${ENTRY}.6tr.fasta
# compress files:
gzip 6tr/${ENTRY}.6tr.fasta
done
```

#### Frame-selection of longest uninterrupted reading frame

```
# List of filenames for the 6-frame translations (.6tr.fasta) where the seq_type == nt:
nt_6tr_path="${MARFERRET_DIR}/scratch/6tr_files.txt"
# Run the 'keep longest frame' script, which outputs the longest reading frame (uninterrupted by stop codons) among the six translations (minimum length of one). 
# If there is a tie for longest coding frame, they will both be kept.
for FASTA in $(cat $nt_6tr_path); do
echo ${FASTA}
keep_longest_frame.py3 -l 1 ${FASTA}
done
```
This script will output a '.6tr.bf1.fasta' file for the 71 translated nucleotide entries. The number of amino acid sequences will be slightly larger than the number of input nucleotide sequences because frames with equal longest-length were output together. 

