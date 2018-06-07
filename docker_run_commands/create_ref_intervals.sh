#!/bin/bash

## Create Reference Intervals for PHG
##
## 


DATA_DIR="ref_intervals"
REF_FA=
GFF=

mkdir -p $DATA_DIR

docker run --rm \
    -v /${DATA_DIR}/:/tempFileDir/data \
    maizegenetics/phg \
    /CreateReferenceIntervals.sh -f -a -e 0