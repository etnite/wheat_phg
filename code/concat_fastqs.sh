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
## *NAME*R1.fastq.gz for forward reads and *NAME*R2.fastq.gz for reverse reads.
##
## The script will match any fastq files matching these patterns, and concatenate
## them into a single fastq file. It doesn't matter if the files are gzipped or
## not.
################################################################################


#### User-defined variables ####

## No trailing slashes!
in_dir="/project/genolabswheatphg/filt_fastqs/wheatCAP_parents"
out_dir="/project/genolabswheatphg/filt_fastqs/wheatCAP_pars_merged"


#### Executable ####

name=$1

## Enable recursive globbing
shopt -s globstar nullglob

## Concatenate forward reads
cat "${in_dir}"/**/*"${name}"*R1.fastq.gz > "${out_dir}"/"${name}"_R1.fastq.gz

## Concatenate reverse reads
cat "${in_dir}"/**/*"${name}"*R2.fastq.gz > "${out_dir}"/"${name}"_R2.fastq.gz


date
