#!/bin/bash

## Create Reference Intervals for PHG
##
## I don't know if this is the behavior of Docker, or of the script that is
## called inside docker here (CreateReferenceIntervals.sh), but it appears that
## both the reference fasta and .gff3 file need to be in the same directory,
## and that the Docker image must be mounted to this directory, and not, for
## instance, one level up in the directory tree.
##
## The three user-defined constants are:
##
##   1) WORK_DIR - directory housing the reference genome fasta and 
##                 corresponding .gff3 annotation file
##   2) REF_FA - Name of reference fasta file
##   3) GFF - Name of .gff3 annotation file
##
## Will probably change this in future to use positional input params.

PROJ_DIR="/workdir/bpward2/wheat_prac_hap_graph"
REF_FA="161010_Chinese_Spring_v1.0_pseudomolecules_parts_renamed.fasta"
GFF="geneModel_v1.gff3"

## Directory containing ref. genome should just be ref/
ref_dir="${PROJ_DIR}"/ref

docker1 run --rm \
    -v "${ref_dir}"/:/tempFileDir/data \
    maizegenetics/phg \
    /CreateReferenceIntervals.sh -f "$REF_FA" -a "$GFF" -e 0

docker1 claim

mv "${ref_dir}"/genomic_intervals* "$PROJ_DIR"/

exit 0;