#!/bin/bash

## This is a simple script to set up the directory structure and data for
## constructing a haplotype graph. This script is just intended for convenience
## and it is assumed that the sources of the input data structures will change,
## requiring periodic modification.

base_dir=$1

mkdir -p base_dir

cd $base_dir

mkdir ref 

## Download v1 RefSeq and .gff3 annotation file
wget -P ref/ http://people.beocat.ksu.edu/~kwjordan/161010_Chinese_Spring_v1.0_pseudomolecules_parts_renamed.fasta
wget -P ref/ http://129.130.89.51/~kjordan/geneModel_v1.gff3