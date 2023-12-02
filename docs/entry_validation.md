The final collection of proteins in each MarFERReT version follows technical validation of each candidate entry for five metrics, and entries failing threshold criteria for these metrics are excluded from the final build (as noted in the 'accepted' column of the metadata.csv and curation.csv files).

Three metrics can be applied to decide on inclusion of new entries:
1. Total sequence count; flag entries with less than 1,200 total sequences
2. Total number of annotated Pfam domains; flag entries with less than 500 Pfam IDs
3. Estimated cross-contamination with ribosomal proteins; flag entries with over 50% sequence cross-contamination
   
The two other metrics were used to decide on entry inclusion for MMETSP re-assembly entries in MarFERReT v1.0 and v1.1  and were added from outside references; these metrics are not necessary for new non-MMETSP candidate entries:
4. MMETSP cross-contamination; flag entries with over 50% sequence cross-contamination reported by [Vlierberghe et al. 2021](https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-021-05717-2) (MMETSP re-assemblies only)
5. MMETSP Ciliate cross-contamination; flag entries estimated cross-contamination reported by [Lasek-Nesselquist and Johnson 2019](https://academic.oup.com/gbe/article/11/11/3218/5610072) (MMETSP re-assemblies only)

### Criteria 1: Total sequence count

Simple metric of the total number of sequences in the candidate entry source .FASTA file. This can be determined in a FASTA-format file with this simple shell command or other tools:

`grep -c ">" entry_sequences.fasta`

Entries with less than 1,200 total sequences are flagged, as noted for the curation.csv file.

### Criteria 2: Total number of annotated Pfam domains

Raw annotation data output by HMMER3 lists can list multiple Pfam domain matches for each query sequence. We select the best-scoring Pfam domain for each sequence and use this annotation for downstream processes and reporting the total number of Pfam annotations. From the raw Pfam annotations (as raw HMMer output for individual entries as .Pfam34.domtblout.tab files or tarballed together in MarFERReT.candidate_entry_Pfam_annotations.tar.gz), use this custom script to select the best Pfam annotation for each sequence from HMMER3 output files:

[best_pfam.py](https://github.com/armbrustlab/marferret/blob/main/scripts/python/best_pfam.py)

The total number of Pfam domains is the size of set of unique Pfam IDs across all Pfam IDs in a candidate entry. Entries with less than 500 unique Pfam domain IDs are flagged for exclusion.

### Criteria 3: Estimated cross-contamination with ribosomal proteins

This is a more complicated step that assesses entries for possible cross-contamination issues. The approach utilizes taxonomic annotation of ribosomal sequences matched to one of 63 Pfam protein domain IDs specific to ribosomal proteins and present in over 90% of candidate entries (referred to in this manuscript as ‘RP63’). To summarize the sub-steps detailed below:

3a. Download reference sequences for ribosomal proteins
3b. Process reference sequences and build a DIAMOND database
3c. Run MarFERReT ribosomal proteins vs reference database
3d. Calculate cross-contamination estimates from LCA annotations

##### 3a. Download sequence references for ribosomal proteins

MarFERReT was designed with the primary purpose of marine metatranscriptome taxonomic annotation, frequently paired with functional annotation through hidden Markov model searches against the known protein family domains in the [Pfam database](https://www.ebi.ac.uk/interpro/entry/pfam/) or similar resources. Criteria 3 is a more complicated analysis to estimate the degree of potential sequence cross-contamination in ribosomal proteins. Entries with low sequence or Pfam count flags from Criteria 1 and 2 above were not considered for this analysis. 

The approach utilizes taxonomic annotation of ribosomal sequences matched to one of 63 Pfam protein domain IDs specific to ribosomal proteins and present in over 90% of candidate entries (referred to in this manuscript as ‘RP63’). The full set of UniProtKB reference sequences associated with each Pfam domain was downloaded from the source through the InterPro database [https://www.ebi.ac.uk/interpro/](https://www.ebi.ac.uk/interpro/). The full set of UniProtKB reference sequences associated with each Pfam domain was downloaded using this script. 

Iterate through the list of Pfam IDs in a call to a custom python script in the bash shell:
``` shell

# **NOTE** this API call is rate-limited and retrieves tens of thousands of # sequences for each Pfam domain. This can take a long time. 
# This can be expedited by registering for an API key with # https://www.ebi.ac.uk

# This loop iterates through each of the 63 Pfam IDs, and passes the 
# $pfam_id to API_Pfam_seq_fasta_download.py

# list of 63 core Pfam domain IDs:
RP63_PFAM_LIST="${MARFERRET_DIR}/data/pfam/ribosomal_proteins/rp63_core_pfams.txt"

# API_Pfam_seq_fasta_download.py is a modified version of the 
# InterPro API call function custom-written to output
# a fasta file with UnitProtKB Pfam sequences identified with this Pfam
# domain (.all_proteins.fasta.gz) and the NCBI taxID to sequence accession file (.all_proteins.mappings.tsv) for DIAMOND database construction.
cd ${MARFERRET_DIR}/data/pfam/ribosomal_proteins/pfam_proteins

for pfam_id in $(cat $RP63_PFAM_LIST); do
	echo "Processing ${pfam_id}" 
	../scripts/python/API_Pfam_seq_fasta_download.py --pfam_id ${pfam_id}
	echo "Sequence count:"
	wc -l ${pfam_id}.reviewed_proteins.mappings.tsv
	echo "~~~~~~~~~~~~~~~"
done

```

##### 3b. Process reference sequences and build a DIAMOND database

After downloading, the RP63 sequences are processed to ensure alignment with the NCBI Taxonomy architecture, to reduce redundancy, and to prepare all files for the correct format for input into the next step with DIAMOND makedb. These processing steps are done in the command shell and documented in this file here:

[RP63_sequence_processing.md](https://github.com/armbrustlab/marferret/blob/main/docs/RP63_sequence_processing.md)

The shell commands detailed in the linked document output these files necessary to create the DIAMOND database:
`Pfam.ribosomal_all_proteins.uniq.faa.gz`
`Pfam.ribosomal_all_proteins.taxonomies.tab.gz`

Create a DIAMOND reference database for the DIAMOND sequence alignment software using the ‘diamond makedb’ command with default parameters and the NCBI database:

``` shell
# generate the diamond library from here:
cd ${MARFERRET_DIR}/data

# retained values from build_diamond_db.sh:
MARFERRET_DIR=$( realpath ../ )
DATA_DIR=${MARFERRET_DIR}/data
TAX_DIR=${DATA_DIR}/diamond/ncbi

TAXONNODES=${TAX_DIR}/nodes.dmp
TAXONNAMES=${TAX_DIR}/names.dmp
VERSION=v1.1

# these values are different than used in build_diamond_db.sh:
REF_PROTEINS="${DATA_DIR}/pfam/ribosomal_proteins/pfam_proteins/Pfam.ribosomal_all_proteins.uniq.faa.gz"
UID2TAXID="${DATA_DIR}/pfam/ribosomal_proteins/pfam_proteins/Pfam.ribosomal_all_proteins.taxonomies.tab.gz"

# output db name:
DMND_DB="${DATA_DIR}/diamond/Pfam.ribosomal_all_proteins.dmnd"

# singularity container is used in this example:
CONTAINER="singularity"
CONTAINER_DIR="${MARFERRET_DIR}/containers"

singularity exec --no-home --bind ${DATA_DIR} \
"${CONTAINER_DIR}/diamond.sif" diamond makedb \
--in ${REF_PROTEINS} --db ${DMND_DB} \
--taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} \
--taxonmap ${UID2TAXID}

# "Database sequences                   6565692
# Database letters                     2179019835
# Accessions in database               6565692
# Entries in accession to taxid file   6565692
# Database accessions mapped to taxid  6565692
# Database sequences mapped to taxid   6565692
# Database hash                        6faaa837ad2db20045ff0f6f5c034a13
# Total time                           32.165000s
# "

```

##### 3c. Run MarFERReT ribosomal proteins vs reference database
``` shell

# we can iterate thru the marferret sequences one at a time
cd ${MARFERRET_DIR}/data/pfam/ribosomal_proteins/mft_ribo_sequences

mkdir lca_full
mkdir lca_full/logs
DMND_DB="${DATA_DIR}/diamond/Pfam.ribosomal_all_proteins.dmnd"
NCORES=4 # this parameter should be adjusted to the cores available
EVALUE="1e-5"

function ribo_diamond {
INPUT_FASTA="${DATA_DIR}/pfam/ribosomal_proteins/mft_ribo_sequences/${ENTRY}.fasta"
LCA_TAB="${DATA_DIR}/pfam/ribosomal_proteins/mft_ribo_sequences/lca_full/${ENTRY}.pfam_full.lca.tab"
echo "Beginning ${ENTRY}"
diamond blastp --tmpdir /scratch/diamond/ -b 100 -c 1 -p $NCORES -d $DMND_DB -e $EVALUE --top 10 -f 102 -q ${INPUT_FASTA} -o ${LCA_TAB} >> lca_full/logs/"${ENTRY}.vs_ribo_ref_full.log" 2>&1
}

for ENTRY in $(ls *fasta | sed 's/.fasta//'); do
echo ${ENTRY}
ribo_diamond ${ENTRY}
done


```

##### 3d. Calculate cross-contamination estimates from LCA annotations

link to:
MarFERReT.RP_validation

Inputs: 	
MAIN INPUT: `.lca.tab`
LIST OF ENTRIES: ex: `mft_ribo_sequences.handles.txt`
METADATA: `taxa.csv`, `metadata.csv`

OUTPUT:  metadata.csv format with QC scores

``` shell

input the hmm file (seq to pfam) 
hmmfile_path = "../../headers/" + this_entry + "_ribo_fasta_headers.csv"

and lca file (seq to taxid)
lcafile_path = this_entry + "_ribosomal.pfam_full.lca.tab"

input metadata.csv
metadat = pd.read_csv("/scratch/marferret/marferret-main/data/MarFERReT.v1.entry_curation.csv") or metadata.csv

input list of entry handles:
with open("mft_ribo_sequences.handles.txt") as f:
	entry_handles = f.read().splitlines()

input taxa.csv
taxa_csv = pd.read_csv("../../Pfam.ribosomal_all_proteins.taxa.csv")
taxa_csv.shape # (112816, 146)

the smaller mft taxa.csv too
mft_taxa_csv = pd.read_csv("/scratch/marferret/marferret-main/data/MarFERReT.v1.taxa.csv")
mft_taxa_csv.shape # (1610, 68)

run the python script:
```


### Visual clustering approach to entry validation

This is not currently used as a metric to decide on acceptance into the quality-controlled MarFERReT build, but we find that hierarchical clustering analysis of Pfam protein-coding content produces an interesting tool for visual identification of aberrant entries. 

The Pfam annotations were used to validate sequence data for candidate entries through hierarchical clustering analysis. We use a custom R script to generate a binary presence/absence matrix of genetic features (Pfam IDs) for each of the entries. The analysis and visualization code can be found in [entry_validation.hclust.R](https://github.com/armbrustlab/marferret/blob/main/scripts/R/entry_validation.hclust.R). 

In the above R code,a binary distance matrix is computed from the presence/absence matrix with the dist() function,and complete hierarchical clustering is conducted with the hclust() function. The ggplot2 and ggtree libraries are used for visualization. From this output, entries are flagged if they did not cluster with their taxonomic lineage or showed other placements indicative of a protein coding inventory incongruent with their known ecological function and taxonomic relationships.
