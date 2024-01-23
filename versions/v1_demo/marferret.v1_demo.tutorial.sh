# MarFERReT demo build tutorial 

# AUTHOR: Ryan Groussman

##################################################
### Demonstration of MarFERReT build processes ###
###  with a small mock data set                ###
##################################################

# This is a streamlined version of the build process designed to test proper processing and production of output files.

# See the main MarFERReT README file for a more in-depth explanation of build process with links to associated files and data sets:
# https://github.com/armbrustlab/marferret/blob/main/README.md#1-cloning-the-marferret-repository


######################################
## 1) Cloning the MarFERReT repository
# Navigate to the demo directory and create your own subdirectory

# This is one way to clone the repository:
gh repo clone armbrustlab/marferret
unzip marferret-main.zip

# This will be the working folder for MarFERReT build processes:
MARFERRET_DIR="${PWD}/marferret-main"

# Let's navigate there
cd $MARFERRET_DIR


######################################
## 2) Collecting and organizing inputs

# 2a. metadata.csv
# To test with the small set; move the original metadata file
# to a separate directory; there hould only be one *metadata.csv file 
# in $MARFERRET_DIR/data/ at any time when running assemble_marferret.sh
mkdir $MARFERRET_DIR/versions/v1.1.1
mv $MARFERRET_DIR/data/MarFERReT.v1.1.1.metadata.csv versions/v1.1.1

# Move the demo metadata.csv file from the /versions/v1_demo/ directory:
mv $MARFERRET_DIR/versions/v1_demo/MarFERReT.v1_demo.metadata.csv .

# 2b. Download source sequences to /source_seqs/

# Retrieve the small mock dataset for the demo
# Sample of 7 entries with 10 sequences each
cd $MARFERRET_DIR/data/source_seqs/
cp $MARFERRET_DIR/versions/v1_demo/demo_mini_entries.tar.gz .
# unpack it:
tar xvfz demo_mini_entries.tar.gz


######################################
## 3) Building software containers

# Complete container setup instructions are available
# in documentation; for the demo we will use Singularity
cd $MARFERRET_DIR/containers/

# Run the Singularity container script
./build_singularity_images.sh


## 4) Running MarFERReT database construction pipeline

# Primary sequence ingestion, standardizing, clustering 
# This script will run within a few minutes 
cd ${MARFERRET_DIR}/scripts
time ./assemble_marferret.sh


#############################################
## 5) Annotating MarFERReT database sequences


## Pfam annotation
# make a pfam directory under the data directory (../data/pfam) 
mkdir ${MARFERRET_DIR}/data/pfam; cd ${MARFERRET_DIR}/data/pfam

# Download the Pfam hmm profiles from their source:
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
# unzip for HMMER3:
gunzip Pfam-A.hmm.gz 

# navigate to the scripts directory and run pfam_annotate.sh 
# (About 1 minute using the small test data)
cd ${MARFERRET_DIR}/scripts
./pfam_annotate.sh 


## Taxonomic annotation with DIAMOND protein-alignment
mkdir ${MARFERRET_DIR}/data/diamond/; mkdir ${MARFERRET_DIR}/data/diamond/ncbi

# the script build_diamond_db.sh will automatically download
# and unzip NCBI taxonomy files if they do not exist 
# Construct a binary DIAMOND database from the MarFERReT build
cd ${MARFERRET_DIR}/scripts
./build_diamond_db.sh

