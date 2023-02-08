A visual diagram of the Part 1 workflow can be found here:
![Part 1 diagram](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/images/diagram1_web.png)

#### Curation of sequence metadata
Manual curation of each entry is necessary to ensure that each sequence entry is standardized with an organism name and an associated NCBI taxonomy ID (tax_id) at the most accurate level possible, if not provided with the source material. (See Methods: Curation of sequence metadata in primary publication for more details). Primary metadata including the MarFERReT entry ID, organism names, NCBI taxID, data type (genome, transcriptome, etc), data source and publication/availability year for current MarFERReT entries are listed here:
[MarFERReT.entry_metadata.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_metadata.v1.csv)

More metadata about the entries including the source organism name, reference publication, original source URL and filename are listed here:
[MarFERReT.entry_source_data.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_source_data.v1.csv)

Optional file listing filenames and relative paths of original and intermediate entry files (for use in building the MarFERReT protein library):
[MarFERReT.entry_paths.v1.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/data/MarFERReT.entry_paths.v1.csv)

#### Gathering
All of the component entries of MarFERReT were aggregated from public and accessible sequence data. Initial sequence data gathering was performed using both command line operations where possible, and through manual downloading and/or naming through a web-based client where otherwise necessary. Command-line acquisition code and instructions for manual downloads are described here: 
[download_source_sequences.md](blob/main/scripts/download_source_sequences.md)

#### Six-frame translation of nucleotide sequences and frame-selection

Sequences gathered in nucleotide alphabet were six-frame translated the protein alphabet, and the longest reading frame (longest uninterrupted stretch of amino acids) was kept for downstream analysis. Instructions for six-frame translation of the 71 nucleotide entries using transeq are here:
[translation_frame_selection.md](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/scripts/translation_frame_selection.md)

#### Functional annotation of protein sequences
Code for downloading Pfam 34.0 and conducting functional annotation for all MarFERReT protein sequences with hmmsearch:
[pfam_annotation.md](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/scripts/pfam_annotation.md)

#### Script for parsing hmmsearch output and selecting the best Pfam for each annotted sequence:
`scripts/best_kofam.ipynb`
[WORK IN PROGRESS]

#### Best-scoring Pfam annotations to MarFERReT on Zenodo:
[LINK TO PFAMS ON ZENODO]
	"MarFERReT_best_pfam.csv"

#### Identification and analysis of Core Transcribed Genes
Core Transcribed Genes (CTGs) were identified from Pfam annotations against MarFERReT eukaryotic transcriptomes using the R programming language. The script for for identifying CTGs for eukaryotic transcriptomes and subsets of major lineages and generating the primary output table and accessory figures can be found here:
[identify_core_transcribed_genes.R](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/scripts/identify_core_transcribed_genes.R)

#### TaxID-level protein clustering
To the reduce sequence redundancy from multiple sequence entries for a single organism, the protein sequences for NCBI taxIDs with more than one entry (source FASTA) were combined and clustered at the 99% amino acid sequence identity threshold. TaxIDs with only a single entry are not clustered. 
[dedupe_and_clustering.md](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/scripts/dedupe_and_clustering.md)

Following this step, the protein sequences from clustered multi-entry taxIDs and the single-entry taxID are combined into a single MarFERRet protein file for subsequent analyses.
