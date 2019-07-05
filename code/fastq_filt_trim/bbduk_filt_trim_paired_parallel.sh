#!/bin/bash

## Quality filtering and adapter trimming using BBDuk
##
## This version of the script runs BBDuk on a single pair of mated fastq files.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for corresponding fastq files using the pattern:
## input_dir/samplename_R[12].fastq.gz
##
## This script is intended to be used with parallelize.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##   1) bbduk parameters are hard coded in the bbduk.sh call below.
##   2) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined constants ####

in_dir="/project/genolabswheatphg/test/BPW_pipeline_test/concat_raw_fastq"
adapt_fasta="/home/brian.ward/repos/wheat_phg/TruSeq_paired_adapters.fa"
out_dir="/project/genolabswheatphg/test/BPW_pipeline_test/filt_fastq"


#### Executable  ####

module load bbtools

echo
echo "Start bbduk_filt_trim_paired_parallel.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"
samp=$1

## Get length of adapters
ad_len=$(head -n 2 "${adapt_fasta}" | tail -n -1 | wc -c)

## Set forward and reverse read fastq files
fq1="${in_dir}"/"${samp}"_R1.fastq.gz
fq2="${in_dir}"/"${samp}"_R2.fastq.gz

## Run BBDuk
## Setting maq will remove entire reads based on their average qual
## maq=10 average error is 10%
## maq=13 ~5%
## maq=15 ~3%
## maq=20 1%
bbduk.sh -Xmx10g in1="${fq1}" in2="${fq2}" \
    out1="${out_dir}"/"${samp}"_R1.fastq.gz out2="${out_dir}"/"${samp}"_R2.fastq.gz \
    ref="${adapt_fasta}" ktrim=r k=$ad_len mink=10 hdist=3 hdist2=1 ftm=5 maq=13 minlen=75 tpe tbo

echo
echo "End time:"
date
