#!/bin/bash

## Index a reference genome for alignment with Hisat2
##
## This script will take a single positional argument - the path to a reference
## genome fasta file, and output a Hisat2 index with the same base name 
## (i.e. the index file names will have the same format as the input reference 
## file name, without the ".fasta" or ".fa" extension).
##
## NOTE: The Hisat2 indexer appears to be single-threaded, unlike the Bowtie2
## indexer
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --job-name="bowtie2-index" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=6 #number of cores/tasks
#SBATCH --time=10:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error

module load hisat2


#### User-defined constants ####

ref_file=$1


#### Executable ####

date

ind_name="${ref_file%.*}"
hisat2-build "${ref_file}" "${ind_name}"

date
