#!/bin/bash -x
# Sample script to run findPath based on fastq from skim sequencing
# and then export a vcf file with imputed variants
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

# Index Pangenome
docker1 run --name cbsu_phg_container_consensus --rm \
	-w / \
	-v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
 	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
	-v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/fastq/ \
	-v ${PHG_DIR}phgSmallSeq/data/PangenomeFasta/:/tempFileDir/outputDir/pangenome/ \
	-v ${PHG_DIR}phgSmallSeq/phgSmallSeq.db:/tempFileDir/outputDir/phgSmallSeq.db \
	-v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
	-t maizegenetics/phg:latest \
	/IndexPangenome.sh phgSmallSeqSequence configSQLiteDocker.txt CONSENSUS 4G 15 10

docker1 run --name small_seq_test_container --rm \
	-w / \
	-v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
 	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
	-v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/fastq/ \
	-v ${PHG_DIR}phgSmallSeq/data/PangenomeFasta/:/tempFileDir/outputDir/pangenome/ \
	-v ${PHG_DIR}phgSmallSeq/phgSmallSeq.db:/tempFileDir/outputDir/phgSmallSeq.db \
	-v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
                    -t maizegenetics/phg:latest \
                    /FindPathMinimap2.sh phgSmallSeqSequence configSQLiteDocker.txt \
                    CONSENSUS CONSENSUS,refRegionGroup \
                    HAP_COUNT_METHOD PATH_METHOD false

echo "Run second docker container - ExportPath.sh"

docker1 run --name cbsu_phg_container_consensus --rm \
	-v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
 	-v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
	-v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/fastq/ \
	-v ${PHG_DIR}phgSmallSeq/data/configSQLiteDocker.txt:/tempFileDir/data/configSQLiteDocker.txt \
	-t maizegenetics/phg:latest \
	/ExportPath.sh configSQLiteDocker.txt CONSENSUS output.vcf

