# Six frame translation was conducted on all sequences ingested in nucleotide-format

# Six-frame translation with the Standard Genetic Code was done on every nucleotide fasta file using transeq: 
transeq -auto -sformat pearson -frame 6 -sequence "$1".fasta -outseq 6tr/"$1".6tr.fasta

