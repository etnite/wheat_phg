#!/bin/bash -x
# Sample script to run findPath based on fastq from skim sequencing
# and then export a vcf file with imputed variants
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

docker1 run --name cbsu_phg_container_consensus --rm \
	-v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
 	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
	-v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/fastq/ \
	-v ${PHG_DIR}phgSmallSeq/phgSmallSeq.db:/tempFileDir/outputDir/phgSmallSeq.db \
	-v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
	-t maizegenetics/phg:latest \
	/FindPath.sh phgSmallSeq.db configSQLiteDocker.txt CONSENSUS Ref.fa HAP_COUNT_METHOD PATH_METHOD

echo "Run second docker container - ExportPath.sh"

docker1 run --name cbsu_phg_container_consensus --rm \
	-v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
 	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
	-v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/fastq/ \
	-v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
	-t maizegenetics/phg:latest \
	/ExportPath.sh configSQLiteDocker.txt CONSENSUS output.vcf
