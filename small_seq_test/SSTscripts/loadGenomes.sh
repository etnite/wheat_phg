#!/bin/bash -x
# This script takes the reference sequence and the reference ranges
# created by createSmallSeqDS.sh and initializes the PHG with them
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

docker1 run --name cbsu_phg_container --rm \
	-v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
	-v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/ \
	-v ${PHG_DIR}phgSmallSeq/answer/:/tempFileDir/answer/ \
	-t maizegenetics/phg:latest \
	/LoadGenomeIntervals.sh configSQLiteDocker.txt Ref.fa anchors.bed Ref_Assembly_load_data.txt true

