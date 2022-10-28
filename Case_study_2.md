## Case Study 2:  Estimating completeness of metatranscriptome taxonomic bins

This Case Study provides an example on how to estimate the completeness of environmental transcriptome bins with taxonomic annotation (Case Study 1) and functional annotation with Pfam 34.0 (ref). The example shown here uses 'species-level' annotations (or lower) for enhanced taxonomic specificity. In summary, the taxonomic and functional annotations are aggregated together and the percentage of lineage-specific core transcribed genes (CTGs) is determined for each species-level environmental taxon bin.

This can be broken down into smaller steps:
1. Pfam annotation of environmental sequences
2. Aggregation of taxonomic and functional annotations
3. Calculation of CTG inventory in transcript bins

### 1. Pfam annotation of environmental sequences

Pfam functional annotation of environmental sequences is performed similarly to annotation of MarFERReT proteins. For downloading Pfam 34.0, see (/scripts/pfam_annotation)[https://github.com/armbrustlab/marine_eukaryote_sequence_database/blob/main/scripts/pfam_annotation.sh]

Declare paths to Pfam profiles, environmental sequences, and TAB output:

`PFAM34_PATH=${MARFERRET_DIR}/pfam/Pfam-A.hmm`
`ENV_SEQS="/mnt/home/user/metatranscriptome/env_assemblies.faa"`
`OUTPUT="/mnt/home/user/metatranscriptome/pfam/env_assemblies.domtblout.tab"`

Use Pfam to annotate environmental contigs against Pfam, using the 'trusted cutoff' score encoded in each hmm profile:

`hmmsearch --cut_tc --domtblout $OUTPUT $PFAM34_PATH $ENV_SEQS`

Convert hmmer output format to csv and select the top-scoring Pfam annotation for each contig. Converted to csv format with a custom script:

`hmmer_domtblout2csv.py -t $OUTPUT -o env_assemblies.pfam.csv`

**WIP - STEPHEN**
Select best-scoring Pfam for each:
"EukRefDB.best_kofam.R"

**WIP - RYAN**

### 2. Aggregation of taxonomic and functional annotations

### 3. Calculation of CTG inventory in transcript bins

identify_core_transcribed_genes.R


