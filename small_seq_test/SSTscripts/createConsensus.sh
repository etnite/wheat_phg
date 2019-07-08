#!/bin/bash -x
# sample script to create consensi for small sequence testing\
# the last 2 parameters are method name for haplotypes (must match
# the value from createHaplotypesAlls.sh script) and the method
# name for this consensus run.
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

docker1 run --name cbsu_phg_container_consensus --rm \
	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
    -v ${PHG_DIR}phgSmallSeq/phgSmallSeq.db:/tempFileDir/outputDir/phgSmallSeq.db \
    -v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
	-v ${PHG_DIR}phgSmallSeq/pangenome/:/tempFileDir/data/outputs/mergedVCFs/ \
	-t maizegenetics/phg:latest \
	/CreateConsensi.sh /tempFileDir/data/configSQLiteDocker.txt Ref.fa GATK_PIPELINE CONSENSUS
