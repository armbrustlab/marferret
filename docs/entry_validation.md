## MarFERReT candidate entry validation by Pfam protein-coding content
We assessed 902 candidate entries for quality and potential contamination issues. Comments, decisions, and supporting evidence for all entries are available in the [MarFERReT entry curation metadata table](https://zenodo.org/record/8011714/files/MarFERReT.v1.entry_curation.csv). 

**Functional annotation of protein sequences**
Candidate entry protein sequences were annotated against the full collection of 19,179 protein family Hidden Markov Models in [Pfam 34.0](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/) using HMMER 3.3 (see [pfam_annotate.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/pfam_annotate.sh) for hmmsearch code and containerization). The highest-stringency cutoff score (‘trusted cutoff’) assigned by Pfam to each hmm profile was used as a minimum score threshold. The general form of the hmmsearch call for each entry (`${ENTRY}`) is:

`hmmsearch --cut_tc --domtblout ${ENTRY}.Pfam34.domtblout.tab Pfam-A.hmm ${ENTRY}.faa`

The raw HMMER output for all MarFERReT entries can be found on the data repository in [MarFERReT.candidate_entry_Pfam_annotations.tar.gz](https://zenodo.org/record/8011714/files/MarFERReT.candidate_entry_Pfam_annotations.tar.gz) (Note: 1.2Gb compressed file). 

If the protein received more than one Pfam match, the best scoring Pfam annotation (highest bitscore) is selected using [best_pfam.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/best_pfam.py) and converted to a CSV format for downstream analyses:

`./best_pfam.py ${ENTRY}.Pfam34.domtblout.tab ${ENTRY}.Pfam34.annotations.csv`

Combine the Pfam.34.annotations.csv files from all entries into a single file, introducing a new file for the identifier column using [combine_pfam_annotations.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/combine_pfam_annotations.py):

`./combine_pfam_annotations.py ${ENTRY_IDS} MarFERReT.v1.aa_seqs.Pfam34.annotations.csv`

Compress the large file:

`gzip MarFERReT.v1.aa_seqs.Pfam34.annotations.csv`

**Hierarchical clustering analysis of Pfam protein-coding content**

The Pfam annotations were used to validate sequence data for 902 candidate entries through hierarchical clustering analysis. Three entries were flagged immediately for low sequence number and coding content; for the remaining 899 we generated a binary presence/absence matrix of 12,549 genetic features (Pfam IDs) for each of the 899 entries. The analysis and visualization of results can be found in [entry_validation.hclust.R](https://github.com/armbrustlab/marferret/blob/main/scripts/R/entry_validation.hclust.R). 

In the above R code,a binary distance matrix is computed from the presence/absence matrix with the dist() function,and complete hierarchical clustering is conducted with the hclust() function. The ggplot2 and ggtree libraries are used for visualization. From this output, entries are flagged if they did not cluster with their taxonomic lineage or showed other placements indicative of a protein coding inventory incongruent with their known ecological function and taxonomic relationships.

