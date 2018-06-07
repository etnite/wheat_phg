#!/bin/bash

## This is a simple script to set up the directory structure and data for
## constructing a haplotype graph. This script is just intended for convenience
## and it is assumed that the sources of the input data structures will change,
## requiring periodic modification.
##
## It takes two positional input parameters - one is the intended path of the
## project directory (absolute), and the second is the relative path of the directory
## within the git repository holding sample configuration files (this is likely
## just "config_files")

proj_dir="$1"
config_dir="$2"

config_dir="$(realpath ${config_dir})"

## Make project directory, and subdirectories
mkdir -p "$proj_dir"
cd "$proj_dir"
mkdir ref config

## Copy the contents of the repository sample config files to the project
## config/ folder
cp "${config_dir}"/* config/

## Download v1 RefSeq and .gff3 annotation file
wget -P ref/ http://people.beocat.ksu.edu/~kwjordan/161010_Chinese_Spring_v1.0_pseudomolecules_parts_renamed.fasta
wget -P ref/ http://129.130.89.51/~kjordan/geneModel_v1.gff3

## Download BAM files