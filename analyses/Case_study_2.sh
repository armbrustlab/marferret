#  Case Study 2: Taxonomic annotations with MarFERReT

#### DIAMOND G*PA vs EukRefDB_v2 dmnd db ####

screen -r diamond

cd /mnt/nfs/projects/ryan/NPacAssemblies_2021/diamond/vs_EukRefDB_v2

EUKREFDB_DMND_DB="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/lc95/dmnd/EukRefDB.v2.dmnd"
NCORES=32
EVALUE="1e-5"

function NPac_diamond {
INPUT_FASTA="/mnt/nfs/projects/ryan/NPacAssemblies_2021/pan-assembly/NPac.${STUDY}.bf100.id99.aa.fasta"
LCA_TAB="/mnt/nfs/projects/ryan/NPacAssemblies_2021/diamond/vs_EukRefDB_v2/NPac.${STUDY}.vs_EukRefDB_v2.lca.tab"
time diamond blastp --tmpdir /scratch/diamond/ -b 100 -c 1 -p $NCORES -d $EUKREFDB_DMND_DB -e $EVALUE --top 10 -f 102 -q ${INPUT_FASTA} -o ${LCA_TAB} >> "${STUDY}.vs_EukRefDB_v2.dmnd.log" 2>&1
}

for STUDY in G1PA G2PA G3PA_UW G3PA_diel; do
NPac_diamond
done
