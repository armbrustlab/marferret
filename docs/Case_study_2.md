## Case Study 2:  Estimating completeness of metatranscriptome taxonomic bins

This Case Study provides an example on how to estimate the completeness of environmental transcriptome bins with taxonomic annotation (Case Study 1) and functional annotation with Pfam 34.0 (ref). The example shown here uses 'species-level' annotations (or lower) for enhanced taxonomic specificity. In summary, the taxonomic and functional annotations are aggregated together and the percentage of lineage-specific core transcribed genes (CTGs) is determined for each species-level environmental taxon bin.

This can be broken down into smaller steps:
1. Pfam annotation of environmental sequences
2. Processing and aggregation of taxonomic and functional annotations
3. Calculation of CTG inventory in transcript bins

### 1. Pfam annotation of environmental sequences

Pfam functional annotation of environmental sequences is performed similarly to annotation of MarFERReT proteins. For instructions on downloading Pfam and annotating sequences with HMMER 3.3, see [pfam_annotate.sh](https://github.com/armbrustlab/marferret/blob/main/scripts/pfam_annotate.sh)

Declare paths to Pfam profiles, environmental sequences, and TAB output:

`PFAM34_PATH=${MARFERRET_DIR}/pfam/Pfam-A.hmm`
`ENV_SEQS=${MARFERRET_DIR}/data/example_env_data/env_assemblies.faa`
`OUTPUT=${MARFERRET_DIR}/data/example_env_data/pfam/env_assemblies.domtblout.tab`

Use Pfam to annotate environmental contigs against Pfam, using the 'trusted cutoff' score encoded in each hmm profile:

`hmmsearch --cut_tc --domtblout $OUTPUT $PFAM34_PATH $ENV_SEQS`

### 2. Processing and aggregation of taxonomic and functional annotations

Pass the raw output from HMMER 3.3 to a script that selects the top-scoring Pfam annotation for each contig, aggregates taxonomic and functional annotations into tables linking the NCBI taxID of a taxon bin to the Pfam IDs associated with the taxID. This produces output tables for Step 3, where the coverage of environmental species bins is estimated from the core transcribed genes (CTGs) derived from [Building Core Transcribed Gene catalog](https://github.com/armbrustlab/marferret/tree/main#6-building-core-transcribed-gene-catalog).

[process_env_annotations.py](https://github.com/armbrustlab/marferret/blob/main/scripts/process_env_annotations.py)

### 3. Calculation of CTG inventory in transcript bins

The output from [process_env_annotations.py](https://github.com/armbrustlab/marferret/blob/main/scripts/process_env_annotations.py) can be passed to these python scripts:

[estimate_completeness.py](https://github.com/armbrustlab/marferret/blob/main/scripts/estimate_completeness.py) will calculate the percent of lineage-dependent CTGs found in each taxon across all samples, as a proxy of completeness.

[estimate_completeness.per_station.py](https://github.com/armbrustlab/marferret/blob/main/scripts/estimate_completeness.per_station.py) does the same calculation of %CTGs for each taxon, but does so for each independent sampling site in the example study. 

