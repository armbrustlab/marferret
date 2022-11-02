# MarFERReT
![marferret logo](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/images/marferret_logo2.png)

## (Marine Functional EukaRyotic Reference Taxa)
### An updated, version controlled, and open-source marine microbial eukaryotic sequence library designed for taxonomic annotation of environmental metatranscriptomes.

## Project Description

The Marine Functional Eukaryotic Reproducible Reference Taxa (MarFERReT) is an  open-source marine eukaryote reference protein sequence library that allows for community-supported expansion over time. MarFERReT was constructed with an emphasis on a comprehensive marine microbial eukaryote reference library for the taxonomic annotation of environmental metatranscriptomes. MarFERRet data was used to identify the core transcribed genes of key marine eukaryote lineages to serve as a metric for estimating the completeness of environmental transcript bins follwoing taxonomic annotations.

The Case Studies presented here are practical examples how MarFERReT can be used by itself or in conjunction with other protein sequence libraries to assign taxonomic identity to environmental sequences using the DIAMOND [fast read aligner], and then to assess the completeness of annotated environmental transcript bins. Lastly, MarFERReT provides a framework for future integration of new microbial eukaryote reference sequences into a growing, versionable sequence library built to increase the scientific community's accessibility and potential in understanding of the function of protists in marine ecosystems.


Primary publication reference:
	[CITATION] Groussman, R.D. et al., ...
	[LINK TO PAPER]

This github repository is associated with a Zenodo repository for data storage:
	[LINK TO ZENODO]

## Table of Contents

The contents of this repo are organized into four parts:
- Part 1: [MarFERReT initial construction](https://github.com/armbrustlab/marine_eukaryote_sequence_database#part-1-marferret-initial-construction)
- Part 2: [Installing and running MarFERRet](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/README.md#part-2-installing-and-running-marferret)
- Part 3: [Case studies]()
- Part 4: [Future MarFERReT releases]()

## Part 1: MarFERReT initial construction

This section details the intial construction of the MarFERReT library and accompanying resources for documentation and replication. A list of important outputs from these methods is included; they are all available in Part 2: Installation. If you want to begin using MarFERReT right away, skip to Part 2.
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

## Part 2: Installing and running MarFERRet
Provide a step-by-step description of how to get the development environment set and running, using the ready-made MarFERReT output generated by the instructions in Part 1.

[STEPHEN STUFF]


## Part 3: Case Study scripts

The Case Use studies here are practical examples how MarFERReT can be used by itself or in conjunction with other protein sequence libraries to assign taxonomic identity to environmental sequences using DIAMOND fast protein alignment, and then to assess the completeness of annotated environmental transcript bins.
A visual diagram of the Case Study workflows can be found here:
![Part 1 diagram](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/images/diagram2_web.png)


[Case Study 1:](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/Case_study_1.md) 
This Case Study shows how MarFERReT can be used to annotate unknown environmental sequences using the DIAMOND fast protein-alignment tool (Buchfink et al., 2015). In summary, a DIAMOND-formatted database is created from sequence data and NCBI Taxonomy information, and used to annotate unknown environmental reads.

[Case Study 2:](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/Case_study_2.md)
This Case Study provides an example on how to estimate the completeness of environmental transcriptome bins with taxonomic annotation (Case Study 1) and functional annotation with Pfam 34.0 (ref). The example shown here uses 'genus-level' annotations (or lower) for enhanced taxonomic specificity. In summary, the taxonomic and functional annotations are aggregated together and the percentage of lineage-specific core transcribed genes (CTGs) is determined for each genus-level environmental taxon bin.



## Part 4: Future MarFERReT releases

MarFERReT was designed to be updated as new microbial eukaryote functional reference sequences are publicly released, with releases identified either through literature reviews, the JGI Genomes On Line Database (GOLD) or through user nominations through the ‘Issues’ request function on github.

New sequences included in MarFERReT will fulfill four requirements: 1) the organism is a marine eukaryote and preferably, a protist, 2) the sequences have been quality-controlled and assembled as transcripts derived from transcriptomes and SATs, or as gene models derived from genomes and SAGs, 3) all sequences are publicly available on a stable repository with an accessible URL, and 4) the organisms should have an associated NCBI taxID at the most-specific taxonomic rank possible.

Users will be able to submit requests for future inclusion through this github repository (link).

When submitting an organism request for future MarFERReT versions, please include the following:
1. Full scientific name of the organism (with strain name if possible)
2. An NCBI taxID of the the organism (as specific as possible, e.g. strain-level)
3. A URL to the location of the assembled source data, with additional instructions if necessary
4. Brief justification for why this organism should be included, e.g. "New SAGs from a clade of marine haptophytes".


## Credits

## License
