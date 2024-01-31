The metadata and code structure underlying MarFERReT were designed to enable the addition of simple and straightforward addition of new reference sequence taxa to existing versions while maintaining the rigor of reproducibility and transparency. 

Different use cases can necessitate different combinations of reference sequence material for comprehensive coverage of taxonomic breadth or better resolution of fine-scale strain diversity. Our lab group has generated hundreds of poly-A filtered (eukaryotic) and ribosomal-rRNA depleted (whole community) environmental metatranscriptomes. We chose to specialize MarFERReT for up-to-date and adaptable focus on marine microbial eukaryotes, and have selected reference taxa primarily for this purpose. 

In practice, we often combine MarFERReT with prokaryote-focused reference sequence libraries for extended domain coverage. See the subsection [Combining MarFERReT with other large reference sequence libraries](https://github.com/armbrustlab/marferret/blob/main/docs/combining_marferret_and_other_references.md#combining-marferret-with-other-large-reference-sequence-libraries) below for code and notes for combining MarFERReT with the MARMICRODB v1.0 prokaryote-focused reference library. 

We plan to generate and publish MarFERReT v1.2 and subsequent releases. We recognize that users may want to add or omit additional reference sequence taxa not included in the current or previous releases. Recommendations for inclusion into future releases are accepted by us through the 'Issues' tab in this repository, but users can also take advantage of this code and metadata framework to include more sequence data and generate derivative database versions. 

In essence, the minimum required to do this is:
1. download the FASTA files for transcriptome assemblies or genome gene models
2. fill in key values in new rows in the metadata.csv table
3. run de novo functional annotations with Pfam using hmmsearch
4. conduct cross-contamination estimates with Pfam (optional)
5. run the ./assemble_marferret.sh script with the new sequence data and metadata

These steps are explained in more detail with code in the subection below:
[Adding new reference sequence entries to MarFERReT individually](https://github.com/armbrustlab/marferret/blob/main/docs/combining_marferret_and_other_references.md#adding-new-reference-sequence-entries-to-marferret-individually)


### Combining MarFERReT with other large reference sequence libraries

MarFERReT can be combined with other domain-focused reference sequence libraries or new reference sequence transcriptomes and genomes to expand taxonomic coverage. One way that we use MarFERReT is in conjunction with a filtered version of the prokaryote-focused MARMICRODB library. Both libraries use NCBI Taxonomy identifiers as their primary classification framework, facilitating compatible annotation approaches. After standardizing data formats, the eukaryoted MarFERReT protein sequences and bacterial MARMICRODB sequences are concatenated together for use in downstream processes. 

#### MarFERReT v1.1 + MARMICRODB v1.0 multi-kingdom marine reference protein sequence library

[marmicrodb_processing.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/marferret_marmicrodb/marmicrodb_processing.sh)
Shell script that describes how to download the MARMICRODB v1.0 files and run a number of accessory scripts to process the MARMICRODB for merging with MarFERReT by removing redundant eukaryotes, environmental MAGs, entries without NCBI taxIDs, and incompatible sequences.
1. Download MARMICRODB database from the Zenodo repository
2. Run filter_marmicrodb_entries.R: Filters out MAGs, eukaryotes, and entries without NCBI taxIDs
3. Run process_marmicrodb_fasta.py: Generate a FASTA file and taxonomy table to combine with MarFERReT
3. Run marmicrodb_remove_numeric_seqs.py: Removes a subset of sequences with numerical values in sequence fields

[merge_marferret_marmicrodb.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/marferret_marmicrodb/merge_marferret_marmicrodb.sh)
Combines the FASTA file and taxonomy tab file from above with the equivalent files from the MarFERReT v1.1 eukaryote sequence library; and contains code for creating an indexed binary database of the combined files for annotation using the DIAMOND fast read alignment program. 
- Includes shell commands to generate Gzip-compressed concatenated FASTA and taxonomy table files. The products are available in this Zenodo repository:
[MarFERReT v1.1 + MARMICRODB v1.0 multi-kingdom marine reference protein sequence library](https://zenodo.org/records/10586950)
- Includes commands to run DIAMOND makedb using above files within the Singularity container. 

### Adding new reference sequence entries to MarFERReT individually

Here's an example on how to get started adding new entries individually broken down into five steps:

**1. Download the FASTA files for transcriptome assemblies or genome gene models**

``` shell

# You will want to create a new subdirectory to get started on this. 

# create a working directory branching off of MarFERReT v1.1
mkdir v1.1-NewTaxa
WORKING_DIR="v1.1-NewTaxa"; cd ${WORKING_DIR}

# new entry download log
mkdir docs
touch docs/new_entry_downloads.md

# download marferret v1.1 metadata.csv and curation.csv files from the
# zenodo data repository:
mkdir ${WORKING_DIR}/data/
cd ${WORKING_DIR}/data/

# download marferret v1.1 metadata.csv and curation.csv files
# note that the downloaded file is renamed to clarify version info
wget -O MarFERReT.v1.1.metadata.csv https://zenodo.org/records/10170983/files/MarFERReT.v1.metadata.csv?download=1 
wget -O MarFERReT.v1.1.curation.csv https://zenodo.org/records/10170983/files/MarFERReT.v1.curation.csv?download=1

# collect your new sequences in FASTA format by downloading or from
# local sources into this directory:
mkdir ${WORKING_DIR}/data/source_seqs

```

**2. Fill in key values in new rows in the metadata.csv table**

As you collect new entries, document how they were retrieved in a log file:
`${WORKING_DIR}/docs/new_entry_downloads.md`

You can see examples of what to write for these logs in examples from previous entries here:
[download_source_sequences.md](https://github.com/armbrustlab/marferret/blob/main/docs/download_source_sequences.md)

As you collect entries, expand the metadata.csv table by adding new rows for each entry you collect, represented by a single source sequence file. This can be done in a separate version as you work, ex: `MarFERReT.v1.1-NewTaxa.metadata.csv`

**3. Run de novo functional annotations with Pfam using hmmsearch**

The new entries should be translated into protein sequence space (if not already) and then given new functional annotations in the same way as the previous entries for uniform functional prediction.  

See this file for an example of how we gave Pfam 34.0 functional domain predictions to MarFERReT protein sequences using HMMER3:
[pfam_annotate.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/pfam_annotate.sh)

**4. Conduct cross-contamination estimates with Pfam (optional)**

Using the Pfam annotations, you can repeat the ribosomal protein cross-contamination estimate method used in MarFERReT v1.1 on new sequence data. We wrote this custom python script to calculate potential cross-contamination estimates and output the `MarFERReT.v1.RP63_QC_estimates.csv` report:
[MarFERReT.RP63_validation.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/MarFERReT.RP63_validation.py)
Note that this file is currently hard-coded to run on existing entries for the MarFERReT v1.1 production; check this repository again for planned modifications to the script allowing it to be called from the command line for individual entries.

See the full documentation in [entry_validation.md](https://github.com/armbrustlab/marferret/blob/main/docs/entry_validation.md) for use of this validation script and other quality control metrics, and the decision tree we used to decide on inclusion into the final quality-controlled v1.1 data product.

**5. Run the assemble_marferret.sh script with the new sequence data and metadata**

When the FASTA sequence data for entries are gathered and the accompanying metadata.csv table is completed with rows for each new entry, this shell script can be used or modified to operate on 'accepted' entries, translate nucleotide data to protein sequences and select the longest coding frames (if applicable), cluster protein sequences for pooled entries sharing taxIDs to reduce intra-taxa redundancy, perform a standardization of sequence names, and generate a final FASTA file(s) with accessory data prepared to create a DIAMOND database for use in annotation: [assemble_marferret.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/assemble_marferret.sh)

This can be ran from the `/scripts/` directory using the command `./assemble_marferret.sh`.

Like MarFERReT, the completed metadata tables and final build products can be deposited and published in a data repository with a stable DOI referencing the contents and creation of your database derivative. Please cite MarFERReT if these derivatives are made.




