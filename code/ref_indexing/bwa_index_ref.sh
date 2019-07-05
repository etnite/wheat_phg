#!/bin/bash

## Index a reference genome for alignment with BWA
##
## This script will take a single positional argument - the path to a reference
## genome fasta file, and output a BWA index with the same base name 
## (i.e. the index file names will have the same format as the input reference 
## file name, without the ".fasta" or ".fa" extension).
##
## NOTE: This script assumes that BWA is installed in a conda environment
## named "bwa". The BWA version must be greater than 0.6 in order to 
## index chromosomes longer than 512Mb. However, at this point, you
## would have to go out of your way to find such an old version...
##
## The bwa index command is single-threaded, and can take a long time (perhaps
## around 8 hours, depending on clock speed) to complete for the hexaploid wheat
## genome
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

module load miniconda


#### User-defined constants ####

ref_file=$1


#### Executable ####

echo
echo "Start bwa_index_ref.sh"
echo "Start time:"
date

source activate bwa

bwa index -a bwtsw "${ref_file}"

source deactivate

echo
echo "End time:"
date
