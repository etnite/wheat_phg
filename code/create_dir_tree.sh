#!/bin/bash

## Create directory structure for PHG
## 
## The PHG requires a specific nested directory structure to operate.
## Docker can create nested directories on the fly, but Singularity can't do this.
## Therefore, we need to just create a directory tree beforehand.


#### User-Defined Constants ####

## Paths to the base directory (the directory which will house the sub-directory tree)
## Path to the reference genome fasta (this will be copied into the sub-dir tree)
## Path to config file (this will also be copied into the sub-dir tree)
BASE_DIR="/workdir/$USER/test_gvcf_load"
REF="/workdir/bpward2/wheat_docker/DataFolders/reference"
CONFIG_FILE="/local/workdir/bpward2/wheat_docker/DataFolders/configSQLite.txt"


#### Executable ####

## Create directory structure for output data
mkdir -p "$BASE_DIR"/answer
mkdir -p "$BASE_DIR"/outputDir/fastq_hap_count
mkdir "$BASE_DIR"/outputDir/pangenome
mkdir "$BASE_DIR"/outputDir/hap_count_best_path
mkdir -p "$BASE_DIR"/data/bam/temp
mkdir "$BASE_DIR"/data/fastq
mkdir "$BASE_DIR"/data/reference
mkdir "$BASE_DIR"/data/gvcfs

## Copy required files
## Reference and gVCFs take up a lot of space, so may want to
## consider deleting these afterwards
cp "$REF" "$BASE_DIR"/data/reference/
cp "$CONFIG_FILE" "$BASE_DIR"/data/config.txt
