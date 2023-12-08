### Part A:  Generating core transcribed genes

This Case Study walks through the approach we call core transcribed genes (CTGs), referring to the actively transcribed genes of recognized function that we can expect to see in a completely sequenced transcriptome in most taxa in a given lineage.

The Pfam annotations of MarFERReT protein sequences are used to identify a set of cross-taxa core transcribed genes that serve as corollary of the BUSCO genome completeness metric, oriented towards marine eukaryotic metatranscriptomes. For any given high-level taxonomic lineage, the CTGs are operationally defined here as the set of Pfam families observed in translated transcriptomes of at least 95% of the species within the given lineage. Only validated entries are used for this analysis. The CTG inventories were identified based on the Pfam 34.0 annotation of 7,514,355 proteins translated from the 654 validated transcriptome and SAT entries (the 146 validated genomic and SAG-sourced entries were not included) of MarFERReT v1.1, and a presence-absence matrix of Pfam functions was generated from the functional annotations of 332 taxa. We derive CTG inventories for all eukaryotes as a whole group, and for nine major marine lineages listed in this custom python script:
[derive_core_transcribed_genes.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/derive_core_transcribed_genes.py)

This requires inputs:
metadata.csv file: `MarFERReT.v1.metadata.csv`
Pfam annotation results in CSV file: `MarFERReT.v1.entry_pfam_sums.csv`
NCBI Taxonomy relationships: `MarFERReT.v1.taxa.csv` (see how to make this file here: [create_ncbi_taxa_csv.md](https://github.com/armbrustlab/marferret/blob/main/docs/create_ncbi_taxa_csv.md)

These are used to call the python script:
``` shell

cd ${MARFERRET_DIR}/scripts/python/

METADATA="${MARFERRET_DIR}/data/MarFERReT.v1.metadata.csv"
PFAM_DATA="${MARFERRET_DIR}/data/MarFERReT.v1.entry_pfam_sums.csv"
TAXA_DATA="${MARFERRET_DIR}/data/MarFERReT.v1.taxa.csv"

./derive_core_transcribed_genes.py --metadata $METADATA \
	--taxa_csv $TAXA_DATA --pfam_dat $TAXA_DATA \
	MarFERReT.v1.core_genes.csv

```

This outputs the CTG catalog as a CSV file: `MarFERReT.v1.core_genes.csv`. 


### Part B:  Using CTGs to assess environmental sequences


This can be broken down into smaller steps:
1. Pfam annotation of environmental sequences
2. Processing and aggregation of taxonomic and functional annotations
3. Calculation of CTG inventory in transcript bins

#### B1. Pfam annotation of environmental sequences

Pfam functional annotation of environmental sequences is performed similarly to annotation of MarFERReT proteins. For instructions on downloading Pfam and annotating sequences with HMMER 3.3, see [pfam_annotate.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/pfam_annotate.sh)

Declare paths to Pfam profiles, environmental sequences, and TAB output:

`PFAM34_PATH=${MARFERRET_DIR}/pfam/Pfam-A.hmm`
`ENV_SEQS=${MARFERRET_DIR}/data/example_env_data/env_assemblies.faa`
`OUTPUT=${MARFERRET_DIR}/data/example_env_data/pfam/env_assemblies.domtblout.tab`

Use Pfam to annotate environmental contigs against Pfam, using the 'trusted cutoff' score encoded in each hmm profile:

`hmmsearch --cut_tc --domtblout $OUTPUT $PFAM34_PATH $ENV_SEQS`

#### B2. Processing and aggregation of taxonomic and functional annotations

Pass the raw output from HMMER 3.3 to a script that selects the top-scoring Pfam annotation for each contig, aggregates taxonomic and functional annotations into tables linking the NCBI taxID of a taxon bin to the Pfam IDs associated with the taxID. This produces output tables for Step 3, where the coverage of environmental species bins is estimated from the core transcribed genes (CTGs) derived from [derive_core_transcribed_genes.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/derive_core_transcribed_genes.py)

#### B3. Calculation of CTG inventory in transcript bins

Generate a binary presence/absence table of Pfam IDs detected in the environmental sample for NCBI taxIDs   (see [derive_core_transcribed_genes.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/derive_core_transcribed_genes.py) for examples of code needed to generate these tables). This processed  Pfam output can be used with the [estimate_completeness.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/estimate_completeness.py) script to calculate the percent of lineage-dependent CTGs found in each taxon across all samples, as a proxy of completeness.

Another example of using putting these CTGs to use  does so for each independent sampling site in the example study rather than the entire study as a whole, facilitating per-site estimates: [estimate_completeness.per_station.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/estimate_completeness.per_station.py) 
