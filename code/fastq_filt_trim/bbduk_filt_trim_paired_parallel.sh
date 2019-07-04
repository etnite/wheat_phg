#!/bin/bash

## Quality filtering and adapter trimming using BBDuk
##
## This version of the script runs BBDuk on a single pair of mated fastq files.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for corresponding fastq files using the pattern:
## input_dir/samplename_R[12].fastq[.gz]
##
## This script is intended to be used with parallelize.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTE: bbduk parameters are hard coded in the bbduk.sh call below.
################################################################################


#### User-defined constants ####

in_dir="/project/genolabswheatphg/raw_data/wheatCAP_parents"
adapt_fasta="/project/genolabswheatphg/Truseq_paired_adapters.fa"
out_dir="/project/genolabswheatphg/filt_fastqs/wheatCAP_parents"


#### Executable  ####

module load bbtools

date
mkdir -p "${out_dir}"
samp=$1

## Get length of adapters
ad_len=$(head -n 2 "${adapt_fasta}" | tail -n -1 | wc -c)

## Set forward and reverse read fastq files
fq1="${in_dir}"/"${samp}"_R1.fastq*
fq2="${in_dir}"/"${samp}"_R2.fastq*

## Run BBDuk
## Setting maq will remove entire reads based on their average qual
## maq=10 average error is 10%
## maq=13 ~5%
## maq=15 ~3%
## maq=20 1%
bbduk.sh -Xmx10g in1="${fq1}" in2="${fq2}" \
    out1="${out_dir}"/"${samp}"_R1.fastq.gz out2="${out_dir}"/"${samp}"_R2.fastq.gz \
    ref="${adapt_fasta}" ktrim=r k=$ad_len mink=10 hdist=3 hdist2=1 ftm=5 maq=13 minlen=75 tpe tbo


date
