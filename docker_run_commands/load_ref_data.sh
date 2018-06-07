#!/bin/bash
set -e

## This script assumes the config.txt file lives in the directory mounted to 
##   /tempFileDir/data
## It assumes the reference fasta (in this case, 
##   Zea_mays.AGPv4.dna.toplevelMtPtv3.fa) lives in the directory mounted to 
##   /tempFileDir/data/reference.
## It assumes the genome intervals file (in the example below, 
##   maizeRefAnchor_intervals_bed) lives in the directory mounted to 
##   /tempFileDir/answer.
## It assumes the genomeData file describing the reference (in the example below, 
##   B73Ref_load_data.txt) lives in the directory mounted to /tempFileDir/data
## It assumes your sqlite database lives in the directory mounted to 
##   /tempFileDir/outputDir/  This in only relevant when running an SQLite 
##   database.  This path shows up in the config file, parameter "db".

WORK_DIR="$1"

mkdir -p "${WORK_DIR}"/database

## Attempt to identify reference fasta in the project ref/ directory
ref_fa="$(basename "${WORK_DIR}"/ref/*.fa*)"
int_bed="$(basename "${WORK_DIR}"/reference_intervals/*gene.bed)"

# You must change "/workdir/user/DockerTuningTests/..." to match your own directory paths
docker1 run --name load_phg_container --rm \
    -v "${WORK_DIR}"/database/:/tempFileDir/outputDir/ \
    -v "${WORK_DIR}"/ref/:/tempFileDir/data/reference/ \
    -v "${WORK_DIR}"/config/:/tempFileDir/data/ \
    -v "${WORK_DIR}"/reference_intervals/:/tempFileDir/answer/ \
    -t maizegenetics/phg \
    /LoadGenomeIntervals.sh config.txt "$ref_fa" "$int_bed" load_reference.txt true

docker1 claim

exit 0;