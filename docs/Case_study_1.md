## Case Study 1:  Taxonomic annotation with DIAMOND protein-alignment

This Case Study shows how MarFERReT can be used to annotate unknown environmental sequences using the DIAMOND fast protein-alignment tool (Buchfink et al., 2015)

This can be broken down into smaller steps:
1. Creating a DIAMOND db
2. Annotating environmental sequences

### 1. Creating a DIAMOND db

The first step is generating a diamond-formatted binary database, combing the actual sequence data with taxonomic information. To make use of the lowest-common-ancest (LCA) algorith with DIAMOND, additional files for NCBI Taxonomy are required.

#### Downloading NCBI Taxonomy files

We use the 'taxtastic' tool to collect NCBI Taxonomy files (taxit v0.9.2). It can be installed on the command line using pip:

`pip install --user -U taxtastic`

Alternately, NCBI Taxonomy files can be downloaded directly from the FTP site:
[https://ftp.ncbi.nih.gov/pub/taxonomy/](https://ftp.ncbi.nih.gov/pub/taxonomy/)

With taxtastic, NCBI Taxonomy database files are automatically downloaded using this command in the current directory:
`taxit new_database`
Unzip the file containing NCBI Taxonomy files:
`unzip taxdmp.zip`

Note: MarFERReT uses NCBI Taxonomy downloaded on 12 October, 2022 at 18:00PST. The NCBI Taxonomy architecture releases regular updates over time, and the NCBI taxIDs, their names and internal classifications are subject to change.

#### Use 'makedb' to create the diamond-formatted database

Declare the relative directory of the MarFERReT base directory. Example:

`MARFERRET_DIR="/mnt/home/user/MarFERReT_v1"`

Declare the paths of all input files:

`cd ${MARFERRET_DIR}/data/dmnd/`
Path to MarFERReT protein fasta (from dedupe_and_clustering.md)
`MARFERRET_FASTA="${MARFERRET_DIR}/data/MarFERReT.v1.proteins.faa"`
Mapping sequence ID to taxID (from uniq_id_and_group_by_taxid.py)
`UID2TAXID="${MARFERRET_DIR}/data/EukRefDB.uid2tax.tab.gz"`

NCBI Taxonomy files:
`TAXONNODES="${MARFERRET_DIR}/data/ncbi/nodes.dmp"
TAXONNAMES="${MARFERRET_DIR}/data/ncbi/names.dmp"`

Name of output diamond db:
`MARFERRET_DMND="${MARFERRET_DIR}/data/dmnd/MarFERReT.v1.dmnd"`

Run DIAMOND makedb with these inputs:
`diamond makedb --in $MARFERRET_FASTA --db ${MARFERRET_DMND} --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} --taxonmap ${UID2TAXID}`


### 2. Annotating environmental sequences

MarFERReT is a protein sequence reference, and is best used to annotate environmental contigs that have already been translated. We encourage you to learn all of the DIAMOND parameters and adjust them to your needs accordingly, and present the parameters used for mapping environmental data with MarFERReT in this paper.

Keep in mind that shotgun sequencing technologies can produce very large volumes of short read data (over hundreds of millions of reads), and annotating short reads directly may result in a significant analysis bottleneck. Assembling short reads into longer contigs generally reduces the number of reads around two orders of magnitude. The environmental sequences used as examples in this case study are assembled and translated from eukaryotic (poly-A) selected metatranscriptomes in the North Pacific (see Lambert et al, 2022 for details).

Declare the path of your environmental protein sequences. Example:

`ENV_SEQS="${DATA_DIR}/metatranscriptome/env_assemblies.faa"`

Default e-value used; change to adjust sensitivity:

`EVALUE="1e-5"`

Using the DIAMOND db created above:

`MARFERRET_DMND="${MARFERRET_DIR}/data/dmnd/MarFERReT.v1.dmnd"`

Define output file for TAB-formatted LCA output (example path):

`LCA_OUT=""${DATA_DIR}/metatranscriptome/vs_MarFERReT_v1/env_seqs.vs_MarFERReT_v1.lca.tab"`

Run diamond blastp in LCA mode (-f 102) using matches within 10% (--top 10) of top alignment score. The "-b 100" and "-c 1" parameters were tuned for system performance and may not be suitable for your machine.

`diamond blastp -b 100 -c 1 -d $DMND_DB -e $EVALUE --top 10 -f 102 -q ${ENV_SEQS} -o ${LCA_OUT}`
