#!/bin/bash

# Use these commands to build a docker image from the Dockerfile
# and save the image to a tarball 

# make and save docker image
docker image build --tag marferret .
docker save marferret | gzip > marferret.tar.gz
