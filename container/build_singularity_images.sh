#!/bin/bash

# Use these command to build a singularity images from
# Docker image tarballs


# build singularity image from tarball
singularity build mmseqs2.sif docker-archive://mmseqs2.tar.gz
singularity build diamond.sif docker-archive://diamond.tar.gz
singularity build hmmer.sif docker-archive://hmmer.tar.gz
singularity build emboss.sif docker-archive://emboss.tar.gz
singularity build marferret-py.sif docker-archive://marferret-py.tar.gz
