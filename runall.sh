# commands to build database from raw sequence by calling each script in order
# and pointing to raw sequence directory and supporting files in data/ directory

MOUNT_DIR=$(realpath data)
IMAGE=container/marferret.sif


singularity exec --bind ${MOUNT_DIR} ${IMAGE} hmmer_domtblout2csv.py --options data/directory/input.tab
