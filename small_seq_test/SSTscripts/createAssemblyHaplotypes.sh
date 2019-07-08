#!/bin/bash -x
# This script loads assemblies into the PHG
# This example is a bit different than how you would normally run the assemblies.
# Normally, a user would process all chromosomes for an assembly within a loop.
# Because our smallSeq dataset has only 1 chromosome, but multiple assemblies, we
# are looping on the assembly name, with the chromosome number remaining constant.
# The referenece fasta also remains constant as it contains only 1 chromosome.
# There are some commented out lines that might be used in a more normal case.
USER=bpward2

# Directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/

#chromList=(1 2 3 4 5 6 7 8 9 10)
assemblyList=(LineA1 LineB1)

#for chrom in "${chromList[@]}"
for assembly in "${assemblyList[@]}"
do

#echo "Starting chrom ${chrom} "
echo "Starting assembly ${assembly} "

docker1 run --name phg_assembly_container_${assembly} --rm \
        -v ${PHG_DIR}phgSmallSeq/:/tempFileDir/outputDir/ \
        -v ${PHG_DIR}phgSmallSeq/ref/:/tempFileDir/data/reference/ \
        -v ${PHG_DIR}phgSmallSeq/data/:/tempFileDir/data/ \
        -v ${PHG_DIR}phgSmallSeq/answer/:/tempFileDir/data/assemblyFasta/ \
        -v ${PHG_DIR}phgSmallSeq/align/:/tempFileDir/outputDir/align/ \
        -t maizegenetics/phg:latest \
        /LoadAssemblyAnchors.sh configSQLiteDocker.txt \
                Ref.fa \
                ${assembly}.fa \
                ${assembly}_Assembly \
                1


echo "Finished chrom  ${assembly} "
done
