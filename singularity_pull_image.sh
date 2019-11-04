#!/bin/bash

## Pull the PHG Docker image using Singularity
##
## Note that the PHG image is actually quite large (about 1.9GB as of this writing), 
## and therefore can take a while to pull. This script just downloads/updates 
## the PHG image, and takes one positional parameter - the path to the directory 
## where the image should be placed.
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="call-vars" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=36:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### Executable ####

if [[ $# -eq 0 ]]; then
    echo "Error - One positional paramter required (path to directory to place image)"
    exit 1;
fi

SINGULARITY_PULLFOLDER="$1"

singularity pull docker://maizegenetics/phg
