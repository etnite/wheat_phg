#!/usr/bin/env python3

'''
Process a bedtools intersect left outer join (-loj)

Brian Ward
brian@brianpward.net
https://github.com/etnite

This script processes the file that is output by bedtools intersect when performing
a left outer join (-loj option). For regions in two .bed files (call them A.bed
and B.bed), performing this operation will create a tab-delimited file with
six columns and no header. The first three columns are the chromosome, start and
end positions from A.bed, while the last three columns are the same
values from B.bed. An intersecting interval from B.bed will be
reported for each interval in A.bed. If there is no match, then a null entry will
be created for B.bed (chromosome represented by a period, start and end positions
represented by -1).

This script will merge together intersecting intervals from A.bed and B.bed. If
there is no intersecting interval in B.bed, then it just prints out the
interval from A.bed.

The script reads from stdin and writes to stdout
'''

import fileinput
import sys


## Read from stdin
for line in fileinput.input():
    
    ## Isolate chrom and positions
    chrom = line.split('\t')[0]
    pos = line.split('\t')[1:3] + line.split('\t')[4:]
    pos = [int(i) for i in pos]
    
    ## If third pos is -1, the B .bed file had a NULL interval. Then just remove
    ## the last two positions (the ones from the B file)
    if pos[2] == -1:
        pos = pos[:2]
    
    lowest = min(pos)
    highest = max(pos)

    ## Print to stdout
    print(chrom + '\t' + str(lowest) + '\t' + str(highest) + '\n')