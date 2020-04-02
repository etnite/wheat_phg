#!/usr/bin/env python3

'''
Process a bedtools intersect left outer join (-loj)

This script processes the file that is output by bedtools intersect when performing
a left outer join (-loj option). For regions in two .bed files (call them A.bed
and B.bed), performing this operation will create a tab-delimited file with
six columns and no header. The first three columns are the chromosome, start and
end positions from A.bed, while the last three columns are the same
values from B.bed. Performing that an intersecting interval from B.bed will be
reported for each interval in A.bed. If there is no match, then a null entry will
be created for B.bed (chromosome represented by a period, start and end positions
represented by -1).

This script will merge together intersecting intervals from A.bed and B.bed. If
there is no intersecting interval in B.bed, then it just prints out the
interval from A.bed.

The script reads in a text file created by bedtools intersect -loj as a positional
parameter
'''

import os
import sys


#### User-Defined Constants

loj_file = sys.argv[1]
#loj_file = '/home/brian/Downloads/SRW_refs_genemodel_left_outer_join.txt'


#### Executable ####

out_file = os.path.splitext(loj_file)[0] + "_parsed.txt"

with open(loj_file, 'r') as in_loj, open(out_file, 'w') as out_bed:
    for line in in_loj:
        chrom = line.split('\t')[0]
        pos = line.split('\t')[1:3] + line.split('\t')[4:]
        pos = [int(i) for i in pos]
        if pos[2] == -1:
            pos = pos[:2]
        
        lowest = min(pos)
        highest = max(pos)

        out_bed.write(chrom + '\t' + str(lowest) + '\t' + str(highest) + '\n')