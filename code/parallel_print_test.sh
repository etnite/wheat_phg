#!/bin/bash

## Test script to print out contents of a file using parallelizer.sh

date
echo "input: $1" 
echo "number of nodes: ${SLURM_NNODES}"
echo "number of cores: ${SLURM_NPROCS}"
echo "job ID: ${SLURM_JOBID}"
