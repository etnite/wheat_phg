#!/bin/bash

## Get SRA samples associated with a BioProject
##
## This script accepts two positional input parameters:
##   1) The accession number of an SRA bioproject
##   2) Output path to write list of samples
##
## The outputted list may need some additional formatting. In addition, a test
## with the v1 wheat hapmap (bioproject SRP032974) showed that some samples had
## multiple associated SRA files - I'm not yet sure why this is the case, or
## how these multiples should be handled.
################################################################################


module load edirect

esearch -db sra -query $1 | efetch --format runinfo | cut -d "," -f 1 | grep -v "Run" | sed '/^[[:space:]]*$/d' | sort | uniq > $2
