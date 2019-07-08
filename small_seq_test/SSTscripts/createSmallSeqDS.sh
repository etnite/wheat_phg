#!/bin/bash -x
# This script creates the dataset in full for the Small Seq Test
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

docker1 run --name cbsu_phg_container --rm \
	-v ${PHG_DIR}configSmallSeq.txt/:/tempFileDir/data/configSmallSeq.txt \
	-v $PHG_DIR:/root/temp/\
	-t maizegenetics/phg:latest \
	/CreateSmallDataSet.sh /tempFileDir/data/configSmallSeq.txt

