#!/bin/bash
set -e

## Convert a .gff3 file to a flattened .bed file
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will take a .gff3 file as input, sort it, isolate only lines that
## represent genes (i.e. get rid of subfeatures), and then will "flatten" it.
## The input .gff3 file can be gzipped or uncompressed.
##
## A second positional input parameter specifies the merge distance - i.e. the
## maximum distance between features (genes in this case) that will
## cause them to be merged together.
##
## In this context, flattening means that overlapping genes will be merged into
## single intervals. The output is a .bed file representing these flattened
## intervals. Bedtools is required to perform the flattening ("merging" in
## bedtools documentation)
##
## NOTES: I don't believe that bedtools merge can read from standard input, so
## a temporary file is created and deleted. Also, note that .gff3 files
## have feature start positions in 1-based coordinates, while these are converted
## to 0-based coordinates in the output .bed file. Bedtools seems to recognize
## gff3 files even if they don't have the .gff3 extension
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="gff2bed" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=1:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### Executable ####

module load bedtools

gff_file=$1
dist=$2

## Remove the ".gff3.gz" for the output name
outname="${gff_file%.*}"
outname=$(echo "${outname}" | sed 's/.gff3//')
outdir=$(dirname "${outname}")

tmpfile=$(mktemp -p "${outdir}")

## Remove comment lines; Get only genes from .gff3 and sort 
zgrep -v "^#" "${gff_file}" |
    grep -P "\tgene\t" |
    sort -k1,1 -k4,4n -k5,5n > "${tmpfile}"

bedtools merge -d $dist -i "${tmpfile}" > "${outname}"_flat.bed

rm "${tmpfile}"
