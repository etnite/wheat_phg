#!/bin/bash
set -e

## Filter BAM files created by Bowtie2
##
## This script is designed to work with arrayer.sh - it takes a single command
## line argument, which is an integer. It uses this integer to select a single
## .bam file from in_dir, then filters the .bam file and writes the output to
## out_dir. 
##
## Notes on isolating "uniquely mapping" reads
##
## This is not a straightforward topic, especially because the concept of
## "uniquely mapping" is ill-defined. Bowtie2 reports a MAPQ score (Q) for each
## alignment which is defined by the probability (p) that an alignment does not
## correspond to the read's true point of origin:
##
##   Q = -10 * log10(p)
##   or:
##   p = 10 ^ -(Q/10)
##
## The value of p is calculated based on the number of alternative alignments,
## and how well these other alignments match the reference. Note that in most
## cases bowtie2 DOES NOT perform an exhaustive search of alignments for each
## read. The user can specify an exhaustive search, but this is ill-advised for
## large genomes (I'm guessing anything larger than bacterial genomes). Therefore,
## the mapq value is based on a sample of the larger set of all possible
## alignments for a given read.
##
## A mapq score of 5 translates to a ~32% chance of misalignment
## A mapq score of 10 translates to a 10% chance
## A mapq score of 20 translates to a 1% chance, etc.
################################################################################ 


#### User-Supplied Constants ####

in_dir="/home/gbg_lab_admin/Array_60TB/wheat_exome_capture/ERSGGL_SRW_alignments/lane_bams"
out_dir="/home/gbg_lab_admin/Array_60TB/wheat_exome_capture/ERSGGL_SRW_alignments/excap_GBS_merged_bams"
mq_thresh=20


#### Executable ####

## Read in array index (integer) - get corresponding .bam file name
arr_ind=$1
bam_file=$(ls -1 *.bam | head -n $arr_ind | tail -n 1)

## Filter input BAM file, sort and output
samtools view -h "${in_dir}"/"${bam_file}" \
    -q $mq_thresh \
    -f 2 |
    samtools sort -O BAM - -o "${out_dir}"/"${bam_file}"

