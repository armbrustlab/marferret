
# six-frame translation:

# Sequences gathered in nucleotide-format gene models
# were translated in six-frames:
transeq -auto -sformat pearson -frame 6 -sequence "$1".fasta -outseq 6tr/"$1".6tr.fasta


# Create a unique id for the database:

REF_DAT="/mnt/nfs/projects/ryan/EukRefDB_2021/metadata/EukRefDB_current.aa_paths.csv"
UID2TAXID="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/test/EukRefDB.uid2tax.tab"
UID2DEFLINE="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/test/EukRefDB.uid2def.csv"
EukRefDB.uniq_id_and_mapping.py -t ${UID2TAXID} -c ${UID2DEFLINE} -r ${REF_DAT}

1. Load in REF_DAT init_ref_dat_dictionary(ref_dat)
2. Use the dictionary to iterate through each Genus_species.
2. Outputs:
	a. New fasta file (Genus_species.fasta )with uid# in defline: >uid#### defline
	b. Mapping file 1: uniq_id to tax_id (tab separated, for DIAMOND input)
	c. Mapping file 2: uniq_id to full defline


# Generate a SPECIES_LIST with all of the unique binomial-named entries (Genus_species)
SPECIES_LIST="/mnt/nfs/projects/ryan/EukRefDB_2021/metadata/EukRefDB.Genus_species_list.txt"


SPECIES_LIST="/mnt/nfs/projects/ryan/EukRefDB_2021/metadata/EukRefDB.Genus_species_list.txt"
MMSEQS_DIR="/mnt/nfs/home/rgrous83/bin/mmseqs/bin"
OUTPUT_DIR="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/lc95"
mkdir ${OUTPUT_DIR}
MIN_SEQ_ID=0.95

function run_linclust {
INPUT_FASTA="${SPECIES}.combined.aa.lensort.fasta"
mkdir ${SPECIES}
mkdir ${SPECIES}/${SPECIES}_tmp
cp /mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/${INPUT_FASTA} ${SPECIES}/
$MMSEQS_DIR/mmseqs createdb ${SPECIES}/${INPUT_FASTA} ${SPECIES}/${SPECIES}.db
$MMSEQS_DIR/mmseqs linclust ${SPECIES}/${SPECIES}.db ${SPECIES}/${SPECIES}.clusters.db ${SPECIES}/${SPECIES}_tmp --min-seq-id ${MIN_SEQ_ID}
$MMSEQS_DIR/mmseqs result2repseq ${SPECIES}/${SPECIES}.db ${SPECIES}/${SPECIES}.clusters.db ${SPECIES}/${SPECIES}.clusters.rep
$MMSEQS_DIR/mmseqs result2flat ${SPECIES}/${SPECIES}.db ${SPECIES}/${SPECIES}.db ${SPECIES}/${SPECIES}.clusters.rep ${SPECIES}/${SPECIES}.aa.lc95.fasta --use-fasta-header
cp ${SPECIES}/${SPECIES}.aa.lc95.fasta ${OUTPUT_DIR}
if [[ -e ${OUTPUT_DIR}/${SPECIES}.aa.lc95.fasta ]]; then
rm -rf /scratch/mmseqs/${SPECIES}; fi
}

# iterate through list and run linclust on all:
for SPECIES in $(cat ${SPECIES_LIST}); do
echo ${SPECIES}
run_linclust
done
