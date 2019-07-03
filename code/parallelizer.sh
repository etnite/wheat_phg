#!/bin/bash

## Parallelize script with GNU Parallel
##
## This script is a "universal parallelizer" - it's intended function is to
## enable the parallel execution of another script WHEN DIFFERENT PROCESSES
## CAN BE RUN INDEPENDENTLY, such as when running the same process on different
## samples, as is common in bioinformatics.
##
## This script must read in an array to use as the parallel iterator. This is
## typically a text file containing a list of some kind with one entry per line.
## Each entry in this list is then fed into whatever script is being called as a
## positional argument. GNU Parallel then handles the parallelization, and
## limits the number of threads. Additional parameters often must be set in the
## script that is being called.
##
## NOTES on using on single machine vs. cluster:
##
## The script is currently set up to run on a cluster using SLURM job management.
## It uses the environmental variable $SLURM_NTASKS to set the number of threads
## based on the user-defined SBATCH lines. If running on a single machine, the
## line "module load parallel" should be commented out, and the number of threads
## constant should be set manually.
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --job-name="parallel-test" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=10  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=6 #number of cores/tasks
#SBATCH --time=00:01:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error

module load parallel


#### User-Defined Constants ####

iter_file=""
script="./my_script.sh"


#### Executable ####

mapfile -t iter < $iter_file

date

parallel -j $SLURM_NTASKS --delay 1 --joblog parallel_run.log srun -n1 $script {} ::: "${iter[@]}"

date
