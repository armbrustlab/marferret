# MarFERReT
![marferret logo](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/images/marferret_logo2.png)

## Marine Functional EukaRyotic Reference Taxa
### A version-controlled, open source library of marine microbial eukaryotic protein sequences for the taxonomic annotation of environmental metatranscriptomes

The Marine Functional EukaRyotic Reference Taxa (MarFERReT) is a version-controlled and open source reference sequence library of marine eukarote proteins that allows for community-supported expansion over time. MarFERReT was constructed for the primary purpose of taxonomic annotation of environmental metatranscriptomes. The case studies included in Part B illustrate how MarFERReT can be used on its own or in combination with other reference libraries for taxonomic and functional annotation, and for estimating the completeness of taxonomic bins.

Primary publication reference:
	[CITATION] Groussman, R.D. et al., ...
	[LINK TO PAPER]

This github repository is associated with a Zenodo repository for data storage:
	[LINK TO ZENODO]

## Table of Contents

The contents of this repo are organized into three parts:
* Part A: [Building MarFERReT](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/README.md#part-a-building-marferret)
* Part B: [Using MarFERRet](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/README.md#part-b-using-marferret)
	* [Case Study 1]()
	* [Case Study 2]()
* Part C: [Future MarFERReT releases](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/README.md#part-c-future-marferret-releases)

## Part A: Building MarFERReT

This section details how to build your own copy of MarFERReT starting from source reference sequences and the scripts stored in this repository. If you want to begin using MarFERReT right away, skip to [Part 2](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/README.md#part-b-using-marferret).

If you're still here, that means you're ready to get into the technical details of building your own copy of the MarFERReT data. This work is broken down into five steps:

1. [Cloning the MarFERReT repository](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/README.md#1-cloning-the-marferret-repository)
1. [Collecting and organizing inputs](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits#2-collecting-and-organizing-inputs)
1. [Building software containers](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits#3-building-software-containers)
1. [Running MarFERReT database construction pipeline](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits#4-running-marferret-database-construction-pipeline)
1. [Annotating MarFERReT database sequences](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits#5-annotating-marferret-database-sequences)

### 1) Cloning the MarFERReT repository

The first step is to copy the MarFERReT pipeline code onto the computer where you intend to build the database. This can be done by [cloning](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) this repository into a suitable directory on your machine. 

### 2) Collecting and organizing inputs

Two sets of input files are required to build MarFERReT: 1) the source reference sequences and 2) a corresponding [metadata file](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/data/MarFERReT.v1.metadata.csv). The source reference sequences will need to be collected from their various public locations, and the metadata file will need to be edited to match the reference sequences. 

#### Source reference sequences

Once you have cloned the MarFERReT repository onto your machine, make a new directory called `source_seqs` under the [`data`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits/data) directory. You will deposit all of the fasta files of the source reference sequences into this directory. Detailed instructions for finding and downloading the source reference sequences used to build MarFERReT v1 can be found [in this document](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/docs/download_source_sequences.md). Before running the MarFERReT pipeline, all of these fasta files should be unzipped. 

#### Metadata file

A metadata file entitled [MarFERReT.v1.metadata.csv](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/data/MarFERReT.v1.metadata.csv) contains important information on each of the source reference sequences used to build the MarFERReT database. Every source reference sequence in the `source_seqs` directory should have a corresponding line in the metadata file with at least the following fields properly filled in:
* `ref_id`: a unique MarFERReT identifier for the reference sequence
* `marferret_name`: a human-readable name for the reference sequence (no spaces or special characters)
* `tax_id`: an NCBI taxonomical identifier
* `source_filename`: this should exactly match the name of the fasta file in `source_seqs` (unzipped)
* `seq_type`: the sequence type of the source fasta -- 'nt' for nucleotide and 'aa' for amino acid
* `aa_fasta`: a name for the standardized fasta file (the convention is 'ref_id' + '_' + 'marferret_name')

### 3) Building software containers

The MarFERReT database construction pipeline is entirely containerized, meaning that you do not need to worry about software dependencies. Additionally, MarFERReT supports both Singularity and Docker containerization, depending on user preference. The necessary containers can be built in two steps:
1. Install either [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) or [Docker](https://docs.docker.com/engine/install/) on your machine. 
1. Navigate to the [`containers`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits/containers) directory and run either the [`build_singularity_images.sh`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/containers/build_singularity_images.sh) or [`build_docker_images.sh`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/containers/build_docker_images.sh) script from the command line.

### 4) Running MarFERReT database construction pipeline

Once the input source reference sequences have been collected, metadata has been organized, and the software containers have been built, you are ready to run the MarFERReT database construction pipeline. Navigate to the [`scripts`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits/scripts) directory and run the [`assemble_marferret.sh`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/scripts/assemble_marferret.sh) script from the command line. You will be prompted to enter either `1` or `2` depending on whether you are using Singularity or Docker containerization. 

The pipeline will take several hours to run, depending on your computer system specifications. When it is done you should find the following outputs in the `data` directory:
* `MarFERReT.v1.proteins.faa.gz` -- MarFERReT protein library
* `MarFERReT.v1.taxonomies.tab.gz` -- taxonomy mapping file required as input for building diamond database
* `MarFERReT.v1.proteins_info.tab.gz` -- mapping file connecting each MarFERReT protein to its originating reference sequence
* `/aa_seq` -- directory with translated & standardized amino acid sequences
* `/taxid_grouped` -- directory with amino acid sequences grouped by taxid
* `/clustered` -- directory with amino acid sequences clustered within taxid

### 5) Annotating MarFERReT database sequences

Information on the functions of the proteins included in MarFERReT can be added in by annotating the sequences with one of the many bioinformatic tools available for functional inference. In this repository we have included a script for annotating the database with [Pfam](https://interpro-documentation.readthedocs.io/en/latest/databases.html#pfam) (now included as a part of the InterPro consortium). 

To annotate MarFERRet, you must first download a copy of the Pfam database of HMM profiles. Make a new directory named `pfam` under the [`data`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits/data) directory. Download into this directory the latest version of Pfam from the [Pfam ftp site](`http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/`). 

Once the Pfam HMM database has been downloaded, navigate to the [`scripts`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits/scripts) directory and run the  [`pfam_annotate.sh`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/scripts/pfam_annotate.sh) script from the command line. In addition to the `data/pfam/Pfam-A.hmm` HMM database, this script requires the  `data/MarFERReT.v1.proteins.faa.gz` file as an input. 

The annotation script can take many days to run, as every protein is compared against every HMM profile. Once it has successfully completed, you should find the following outputs in the `data` directory:
* `MarFERReT.v1.best_pfam_annotations.csv.gz` -- a summary of the best Pfam annotation for each MarFERReT reference protein
* `pfam/MarFERReT.${VERSION}.pfam.domtblout.tab.gz` -- the complete set of Pfam annotations against each MarFERReT reference protein

## Part B: Using MarFERReT

The primary intended use of MarFERReT is the taxonomical annotation of marine metatranscriptomic datasets. This can be done without building your own copy of the database. Instead, the MarFERReT v1 database files necessary for annotaiton can be downloaded from [Zenodo]().

One means of performing this taxonomical annotation is to search MarFERReT for the closest matches to your data sequences. There are many bioinformatic tools available for this type of sequence alignment. One such popular tool for high performance sequence alignment of big datasets is [DIAMOND](https://github.com/bbuchfink/diamond). In the [`scripts`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/tree/sbedits/scripts) directory of this repository, we have included a script named [`build_diamond_db.sh`](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/sbedits/scripts/build_diamond_db.sh) for using DIAMOND in combination with MarFERReT. This script requires the `MarFERReT.v1.proteins.faa.gz` and `MarFERReT.v1.taxonomies.tab.gz` as inputs in the `data` directory.

### Case Studies

The case studies presented here are practical examples how MarFERReT can be used by itself or in conjunction with other protein sequence libraries to assign taxonomic identity to environmental sequences using DIAMOND fast protein alignment, and then to assess the completeness of annotated environmental transcript bins.

A visual diagram of the case study analyses can be found here:
![Part 1 diagram](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/images/diagram2_web.png)

[Case Study 1:](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/Case_study_1.md) 
This Case Study shows how MarFERReT can be used to annotate unknown environmental sequences using the DIAMOND fast protein-alignment tool (Buchfink et al., 2015). In summary, a DIAMOND-formatted database is created from sequence data and NCBI Taxonomy information, and used to annotate unknown environmental reads.

[Case Study 2:](https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/Case_study_2.md)
This Case Study provides an example on how to estimate the completeness of environmental transcriptome bins with taxonomic annotation (Case Study 1) and functional annotation with Pfam 34.0 (ref). The example shown here uses 'genus-level' annotations (or lower) for enhanced taxonomic specificity. In summary, the taxonomic and functional annotations are aggregated together and the percentage of lineage-specific core transcribed genes (CTGs) is determined for each genus-level environmental taxon bin.

## Part C: Future MarFERReT releases

MarFERReT was designed to be updated as new microbial eukaryote functional reference sequences are publicly released, with releases identified either through literature reviews, the JGI Genomes On Line Database (GOLD) or through user nominations via the ‘Issues’ request function on github.

New sequences included in MarFERReT must fulfill four requirements: 1) the organism is a marine eukaryote, and preferably a protist, 2) the sequences have been quality-controlled and assembled as transcripts derived from transcriptomes and SATs, or as gene models derived from genomes and SAGs, 3) all sequences are publicly available on a stable repository with an accessible URL, and 4) the organisms should have an associated NCBI taxID at the most-specific taxonomic rank possible.

Users will be able to submit requests for future inclusion through this github repository.

When submitting an organism request for future MarFERReT versions, please include the following:
1. Full scientific name of the organism (with strain name if possible)
1. An NCBI taxID of the the organism (as specific as possible, e.g. strain-level)
1. A URL to the location of the assembled source data, with additional instructions if necessary
1. Brief justification for why this organism should be included, e.g. "New SAGs from a clade of marine haptophytes".

## Credits

## License
