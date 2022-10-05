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

Note: MarFERReT uses NCBI Taxonomy downloaded on 13 January, 2022 at 18:00PST. The NCBI Taxonomy architecture changes over time and NCBI taxIDs, their names and classifications are subject to change.

#### Use 'makedb' to create the diamond-formatted database

Declare the relative directory of the MarFERReT base directory. Example:

`MARFERRET_DIR="/mnt/home/user/MarFERReT_v1"`

Declare the paths of all input files:

```
cd ${MARFERRET_DIR}/dmnd/

# (note, formerly "EukRefDB.combined.lc95.uid.faa")
# MarFERReT protein file:
IN_FASTA="${MARFERRET_DIR}/data/seqs/MarFERReT.dmnd.aa.faa"
# TAB file linking protein sequence to taxID
TAXONMAP="${MARFERRET_DIR}/data/MarFERReT.tax_mapping.tab.gz"
# NCBI Taxonomy files:
TAXONNODES="${MARFERRET_DIR}/NCBI_db/nodes.dmp"
TAXONNAMES="${MARFERRET_DIR}/NCBI_db/names.dmp"
# Name of output diamond db:
DMND_DB="${MARFERRET_DIR}/dmnd/MarFERReT.dmnd"
time diamond makedb --in $IN_FASTA --db ${EUKREFDB_DMND_DB} --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} --taxonmap ${TAXONMAP}
```

### 2. Annotating environmental sequences

MarFERReT is a protein sequence reference, and is best used to annotate environmental contigs that have already been translated. We encourage you to learn all of the DIAMOND parameters and adjust them to your needs accordingly, and present the parameters used for mapping environmental data with MarFERReT in this paper.

Keep in mind that shotgun sequencing technologies can produce very large volumes of short read data (over hundreds of millions of reads), and annotating short reads directly may result in a significant analysis bottleneck. Assembling short reads into longer contigs generally reduces the number of reads around two orders of magnitude. The environmental sequences used as examples in this case study are assembled and translated from eukaryotic (poly-A) selected metatranscriptomes in the North Pacific (see Lambert et al, 2022 for details).

Declare the path of your environmental protein sequences. Example:

`ENV_SEQS="/mnt/home/user/metatranscriptome/env_assemblies.faa"`

Code below is designed to run in the diamond subdirectory:

`cd ${MARFERRET_DIR}/dmnd/`

Default e-value used; change to adjust sensitivity:

`EVALUE="1e-5"`
Using the DIAMOND db created above:

`DMND_DB="${MARFERRET_DIR}/dmnd/MarFERReT.dmnd"`

Output file for LCA analysis (TAB format):
`LCA_OUT="/mnt/nfs/projects/ryan/NPacAssemblies_2021/diamond/vs_EukRefDB_v2/NPac.${STUDY}.vs_EukRefDB_v2.lca.tab"`

Run diamond blastp in LCA mode (-f 102) using matches within 10% (--top 10) of top alignment score. The "-b 100" and "-c 1" parameters were tuned for system performance and may not be suitable for your machine.
`diamond blastp -b 100 -c 1 -d $DMND_DB -e $EVALUE --top 10 -f 102 -q ${ENV_SEQS} -o ${LCA_OUT}`

#### References

BUCHFINK REF 
