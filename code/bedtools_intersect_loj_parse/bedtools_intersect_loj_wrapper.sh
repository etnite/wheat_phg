#!/bin/bash

## Wrapper to parse output of bedtools intersect -loj
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script takes two .bed files as input. It then performs bedtools intersect
## using the -loj (Left outer join) option. This creates a file with six columns:
## the chromosome, start, and end positions from the two input .bed files.
## For features that are in the first .bed file (a_bed_file) that don't overlap
## any feature in the second .bed file (b_bed_file), a null feature is written
## for b_bed_file, where the chromosome is listed as ".", and start and end
## positions are listed as -1. This output is streamed through the associated
## parse_bedtools_intersect_loj.py, which finds the minimum and maximum positions
## for each intersection, and outputs a three-column .bed file.
##
## This .bed file is then run through bedtools merge to consolidate any repeated/
## overlapping intervals. (This happens when, for instance, there are multiple
## features in a_bed_file overlapping one feature in b_bed_file). The script is
## currently configured to also remove the wheat Un "chromosome" of unaligned
## contigs.
################################################################################


#### User-Defined Constants ####

## Paths to the two input .bed file and the output .bed file
a_bed_file="/home/brian/Downloads/A_test.bed"
b_bed_file="/home/brian/Downloads/B_test.bed"
out_bed_file="/home/brian/Downloads/test_loj.bed"

## Maximum merge distance for bedtools merge
## Setting to 0 will only merge overlapping/bookended features
merge_dist=0


#### Executable ####

#module load bedtools

out_dir=$(dirname "$out_bed_file")
mkdir -p "$out_dir"
tmpfile=$(mktemp -p "$out_dir")

## Perform the intersection, parse the output and sort it for good measure
bedtools intersect -sorted -loj -a "$a_bed_file" -b "$b_bed_file" |
	./parse_bedtools_intersect_loj.py |
	sort -k 1,1 -k2,2n > "$tmpfile"

## Unfortunately bedtools merge can't read from stdin
## We are also removing the Un "chromosome" here
bedtools merge -d $merge_dist -i "$tmpfile" |
	grep -v "^Un" > "$out_bed_file"

rm "$tmpfile"

