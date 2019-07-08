#!/bin/bash -x
# This script creates haplotypes on the basis of fastq rather than assemblies
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

taxonList=(Ref LineA LineB RefA1 LineA1 LineB1)

for taxon in "${taxonList[@]}"
do
mkdir -p /workdir/${USER}/PHG_training_tests/dockerOutput/gvcfOut/${taxon}/
mkdir -p /workdir/${USER}/PHG_training_tests/dockerOutput/gvcfOutFilter/${taxon}/

docker1 run --name cbsu_phg_container_${taxon} --rm \
	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
    -v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/fastq/ \
    -v ${PHG_DIR}phgSmallSeq/phgSmallSeq.db:/tempFileDir/outputDir/phgSmallSeq.db \
	-v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
	-v ${PHG_DIR}phgSmallSeq/align/:/tempFileDir/data/outputs/gvcfs/ \
	-t maizegenetics/phg:latest \
	/CreateHaplotypes.sh /tempFileDir/data/configSQLiteDocker.txt \
			  ${taxon} \
			  single \
			  GATK_PIPELINE \
			  ${taxon}_R1.fastq

done
