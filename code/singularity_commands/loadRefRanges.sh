#!/bin/bash

## Load References Ranges using Singularity Image
##
## This script assumes the following: 
##   * config.txt file lives in the directory mounted to /tempFileDir/data
##   * reference fasta lives in the directory mounted to /tempFileDir/data/reference.
##   * the genome intervals file lives in the directory mounted to /tempFileDir/answer.
##   * the genomeData file describing the reference (in the example below, wheat_load_data.txt) lives in the directory mounted to /tempFileDir/data
##   * sqlite database lives in the directory mounted to /tempFileDir/outputDir/  This in only relevant when running an SQLite database.  This path shows up in the config file, parameter "db".
#####################################################################################################


#### User-Defined Constants

BASE_DIR=""

## These parameters require JUST THE BASENAME of each file - the script assumes
## each file is located in the paths described above
CONFIG_FILE="config.txt"
REF=""
INTERVALS_FILE=""
GENOME_FILE=""


#### Executable ####

singularity run \
        -B /workdir/${USER}/wheat_docker/DockerOutput/:/tempFileDir/outputDir/ \
        -B /workdir/${USER}/wheat_docker/DataFolders/:/tempFileDir/data \
        -B /workdir/${USER}/wheat_docker/DataFolders/:/tempFileDir/answer/ \
        /workdir/${USER}/wheat.simg \
        /LoadGenomeIntervals.sh "$CONFIG_FILE" "$REF" "$INTERVALS_FILE" "$GENOME_FILE" true
