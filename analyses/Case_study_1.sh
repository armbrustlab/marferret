# Case Study 1: Building the DIAMOND db


### Creating a DIAMOND db from MarFERReT alone:


# generate the diamond_db file with these inputs:

screen -r refdb
cd /mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/lc95/dmnd/

EUKREFDB_FASTA="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/lc95/EukRefDB.combined.lc95.uid.faa"
TAXONNODES="/mnt/nfs/projects/ryan/EukRefDB_2021/NCBI_db/nodes.dmp"
TAXONNAMES="/mnt/nfs/projects/ryan/EukRefDB_2021/NCBI_db/names.dmp"
TAXONMAP="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/EukRefDB.uid2tax3.tab.gz"
EUKREFDB_DMND_DB="/mnt/nfs/projects/ryan/EukRefDB_2021/combined_seqs/lc95/dmnd/EukRefDB.v2.dmnd"
time diamond makedb --threads 32 --in $EUKREFDB_FASTA --db ${EUKREFDB_DMND_DB} --taxonnodes ${TAXONNODES} --taxonnames ${TAXONNAMES} --taxonmap ${TAXONMAP}


### Joining MarFERReT with another reference database
