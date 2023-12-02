RP63 sequences are processed to ensure alignment with the NCBI Taxonomy architecture, to reduce redundancy, and to prepare all files for the correct format for input into the next step with DIAMOND makedb. 

##### Step 1: Ensure alignment with NCBI Taxonomy

``` shell
# Retrieve unique IDs from the mapping file:
cd ${MARFERRET_DIR}/data/pfam/ribosomal_proteins/pfam_proteins
zcat Pfam.ribosomal_all_proteins.taxonomies.tab.gz | tail -n +2 | awk -F"\t" '{print $3}' | uniq | sort | uniq > Pfam.ribosomal_all_proteins.uniq_taxIDs.txt

# Check the total count:
wc -l Pfam.ribosomal_all_proteins.uniq_taxIDs.txt
#	"79995 Pfam.ribosomal_all_proteins.uniq_taxIDs.txt"
	
# Run these taxIDs through the taxit tool to identify
# IDs that are not compatible with current NCBI taxonomy database
TAX_LIST="Pfam.ribosomal_all_proteins.uniq_taxIDs.txt"
TAX_DB="${MARFERRET_DIR}/data/diamond/ncbi/ncbi_taxonomy.db"
taxit taxtable -f ${TAX_LIST} -o ../Pfam.ribosomal_all_proteins.taxa.csv ${TAX_DB}

# This identifies a handful of PROBLEM TAXIDS that are a small
# subset of the tens of thousands in total and can be removed.
# ['1031566', '1055723', '11103', '1133492', '121791', '1345260', '1407013', '1430529', '1497955', '1715159', '1740162', '1783512', '182327', '183483', '1847729', '1875017', '1919277', '196180', '2024845', '2026160', '2045443', '228257', '23094', '2448762', '2448765', '2583024', '2654173', '2663858', '2675300', '2679994', '2718658', '2740746', '2751884', '2771433', '2843217', '2849021', '2978340', '2995172', '2995173', '342610', '433476', '444923', '458840', '477290', '53027', '59244', '63330', '645273', '645274', '648567', '66091', '66898', '67334', '80870', '866775', '93032', '938176', '977743', '977746', '977749']

# make a list of the 60 IDs to remove in a file:
REMOVE_IDS="Pfam.ribosomal_all_proteins.removed_taxIDs.txt"

# and iterate through the list of taxIDs:
TAXIDS_IN="Pfam.ribosomal_all_proteins.uniq_taxIDs.txt"
grep -Fwv -f ${REMOVE_IDS} ${TAXIDS_IN} > Pfam.ribosomal_all_proteins.uniq_taxIDs2.txt

# to create a list with the 60 problematic IDs removed
wc -l Pfam.ribosomal_all_proteins.uniq_taxIDs2.txt
	"79935 Pfam.ribosomal_all_proteins.uniq_taxIDs2.txt"

# Now try taxit again:
TAX_LIST="Pfam.ribosomal_all_proteins.uniq_taxIDs2.txt"
TAX_DB="${MARFERRET_DIR}/data/diamond/ncbi/ncbi_taxonomy.db"
taxit taxtable -f ${TAX_LIST} -o ../Pfam.ribosomal_all_proteins.taxa.csv ${TAX_DB}

# If the taxIDs are aligned, taxit produces the taxa.csv table for all taxa represented in the Pfam sequences:
TAXA_CSV="${MARFERRET_DIR}/data/pfam/ribosomal_proteins/Pfam.ribosomal_all_proteins.taxa.csv"

```

##### Step 2: Reduce sequence redundancy
``` shell 

cd ${MARFERRET_DIR}/data/pfam/ribosomal_proteins/pfam_proteins

# Reduce sequence redundancy from annotation overlaps
# There are this many protein sequences;
cat *fasta | grep -c ">" # 7068579
# and this many unique protein sequences:
cat *.fasta | grep ">" | sed 's/>//' | sort | uniq | wc # 6565692

# concatenate reference proteins together:
cat *.all_proteins.fasta >> Pfam.ribosomal_all_proteins.faa

# filter to only uniques:
seqmagick convert --deduplicate-taxa Pfam.ribosomal_all_proteins.faa  Pfam.ribosomal_all_proteins.uniq.faa
# The unique sequence count in the fasta file after filtering 
grep -c ">" Pfam.ribosomal_all_proteins.uniq.faa # 6565692
# compress it:
gzip Pfam.ribosomal_all_proteins.uniq.faa
# 941M Oct 13 15:15 Pfam.ribosomal_all_proteins.uniq.faa.gz

# remove/zip intermediate files:
rm Pfam.ribosomal_all_proteins.faa
for fasta in $(ls *.all_proteins.fasta); do
echo $fasta
gzip $fasta
done

# finish constructing the mapping file:
# add this header to the tsv file first
echo "accession	accession.version	taxid	gi" > Pfam.ribosomal_all_proteins.taxonomies.tab
# make sure the correct tab format is preserved

# concatenate the rest, remove redundancy:
cat *.all_proteins.mappings.tsv | sort | uniq >> Pfam.ribosomal_all_proteins.taxonomies.tab
wc -l Pfam.ribosomal_all_proteins.taxonomies.tab # 6565693

# compress the mapping file:
gzip Pfam.ribosomal_all_proteins.taxonomies.tab

```
