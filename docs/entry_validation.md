## Validation of MarFERReT candidate entries

### Necessary input files:
#### entry curation metadata:
"MarFERReT.v1.entry_curation.csv"
#### NCBI taxonomy relationship file:
"MarFERReT.v1.taxa.csv"
#### Pfam annotations for all entries, unclustered
"MarFERReT.candidate_entry_Pfam_annotations.tar.gz "

#### use these wrapper functions/scripts to get from the raw file:
"MarFERReT.candidate_entry_Pfam_annotations.tar.gz"
#### and process the raw data to get a CSV file with the annotations:
"MarFERReT.v1.aa_seqs.Pfam34.annotations.csv.gz"
    
#### build the presence/absence matrix, binary distance matrix, and conduct
#### hierarchical clustering in this R script, along with analytical figures using ggplot2. 
"entry_validation.hclust.R"



