# MarFERReT

## (Marine Functional EukaRyotic Reference Taxa)
### An updated, version controlled, and open-source marine microbial eukaryotic sequence library designed for taxonomic annotation of environmental metatranscriptomes.

## Project Description

The Marine Functional Eukaryotic Reproducible Reference Taxa (MarFERReT) is an  open-source marine eukaryote reference protein sequence library that allows for community-supported expansion over time. MarFERReT was constructed with an emphasis on a comprehensive marine microbial eukaryote reference library for the taxonomic annotation of environmental metatranscriptomes. MarFERRet data was used to identify the core transcribed genes of key marine eukaryote lineages to serve as a metric for estimating the completeness of environmental transcript bins follwoing taxonomic annotations.

The Case Use studies here are practical examples how MarFERReT can be used by itself or in conjunction with other protein sequence libraries to assign taxonomic identity to environmental sequences using the DIAMOND [fast read aligner], and then to assess the completeness of annotated environmental transcript bins. Lastly, MarFERReT provides a framework for future integration of new microbial eukaryote reference sequences into a growing, versionable sequence library built to increase the scientific community's accessibility and potential in understanding of the function of protists in marine ecosystems.


Primary publication reference:
	[CITATION] Groussman, R.D. et al., ...
	[LINK TO PAPER]

This github repository is associated with a Zenodo repository for data storage:
	[LINK TO ZENODO]

## Table of Contents

The contents of this repo are organized into three main categories:
- Part 1: MarFERReT initial construction (optional)
- Part 2: Installing and running MarFERRet (the primary use of MarFERReT)
- Part 3: Case Study scripts (optional analysis scripts)
- Part 4: Future MarFERReT releases

## Part 1: MarFERReT initial construction

This section details the intial construction of the MarFERReT library and accompanying resources for documentation and replication. A list of important outputs from these methods is included; they are all available in Part 2: Installation. If you want to begin using MarFERReT right away, skip to Part 2.
A visual diagram of the Part 1 workflow can be found here:
[LINK TO DIAGRAM]

#### Gathering
All of the component entries of MarFERReT were aggregated from public and accessible sequence data. Initial sequence data gathering was performed using both command line operations where possible, and through manual downloading and/or naming through a web-based client where otherwise necessary. Command-line acquisition code is found here:
`/scripts/download_source_sequences.sh`

Manually-downloaded data is described in as much detail as possible here:
`scripts/webclient_source_sequences.md`

#### Curation of sequence metadata
After downloading the raw source material, manual curation is necessary to ensure that each sequence entry is standardized with an organismal name and an associated NCBI taxonomy ID (tax_id), if not provided with the source material. (See Methods: Curation of sequence metadata in primary publication for more details). The results of this manual curation can be found in the primary MarFERReT entry metadata file:
`MarFERReT_entries.v1.csv` [WIP]

#### Six-frame translation of nucleotide sequences
Code for six-frame translation using EMBOSS is here:
`scripts/translation.sh`

#### Functional annotation of protein sequences
Code for downloading Pfam 34.0 and conducting functional annotation of MarFERReT protein sequences with hmmsearch:
`scripts/pfam_annotation.sh`

#### Script for parsing hmmsearch output and selecting the best Pfam for each annotted sequence:
`scripts/best_kofam.ipynb`
[WORK IN PROGRESS]

#### Best-scoring Pfam annotations to MarFERReT on Zenodo:
[LINK TO PFAMS ON ZENODO]
	"MarFERReT_best_pfam.csv"

#### Identification and analysis of Core Transcribed Genes
Code for deriving Core Transcribed Genes in MarFERReT eukaryotic transcriptomes, conducted in R:
`scripts/identify_core_transcribed_genes.R`

#### Species-level protein clustering
Code for species-level protein clustering
`scripts/clustering.sh`


## Part 2: Installing and running MarFERRet
Provide a step-by-step description of how to get the development environment set and running.

[STEPHEN STUFF]


## Part 3: Case Study scripts

The Case Use studies here are practical examples how MarFERReT can be used by itself or in conjunction with other protein sequence libraries to assign taxonomic identity to environmental sequences using DIAMOND fast protein alignment, and then to assess the completeness of annotated environmental transcript bins.

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
