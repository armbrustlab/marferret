#!/bin/bash

# Use this command to build a singularity image from a 
# Docker image tarball


# build singularity image from tarball
singularity build marferret.sif docker-archive://marferret.tar.gz
