#!/bin/bash

## Concatenate fastq files
##
## This script takes as input a single string (typically a genotype NAME). It
## will typically be used with parallelizer.sh to run in independent, parallel
## fashion.
##
## The user must define an input directory, and an output directory.
##
## The fastq files should be in the input directory, named in the format 
## *NAME*_R1_*.fastq.gz for forward reads and *NAME*_R2_*.fastq.gz for reverse reads.
##
## The script will match any fastq files matching these patterns, and concatenate
## them into a single fastq file. It doesn't matter if the files are gzipped or
## not.
##
## NOTES:
##   1) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
##
################################################################################


#### User-defined variables ####

## No trailing slashes!
in_dir="/project/genolabswheatphg/raw_data/wheatCAP_parents/separated_fq"
out_dir="/project/genolabswheatphg/test/BPW_pipeline_test/concat_raw_fastq"


#### Executable ####

echo
echo "Start concat_fastqs.sh"
echo "Start time:"
date

name=$1
mkdir -p "${out_dir}"

## Enable recursive globbing
shopt -s globstar nullglob

## Concatenate forward reads
cat "${in_dir}"/**/*"${name}"*_R1_*.fastq.gz > "${out_dir}"/"${name}"_R1.fastq.gz

## Concatenate reverse reads
cat "${in_dir}"/**/*"${name}"*_R2_*.fastq.gz > "${out_dir}"/"${name}"_R2.fastq.gz

echo
echo "End time:"
date
