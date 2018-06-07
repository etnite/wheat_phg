#!/bin/bash

mkdir -p ref ref_intervals

## Download v1 RefSeq and .gff3 annotation file
wget -P ref/ http://people.beocat.ksu.edu/~kwjordan/161010_Chinese_Spring_v1.0_pseudomolecules_parts_renamed.fasta
wget -P ref/ http://129.130.89.51/~kjordan/geneModel_v1.gff3