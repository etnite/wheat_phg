#! /bin/bash

## Create PHG haplotypes from gVCF files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This wrapper script will proceed through the steps necessary to create a PHG
## from either read data (fastq files), or alignments (.bam or gVCF files)
##
## The steps performed are:
##   1) LoadGenomeIntervals
##   2) CreateHaplotypes
##
## Loading from gVCF files is the fastest way to load haplotypes into the PHG.
## Typically I use gVCFs as inputs, as the conversion from BAM files to gVCFs
## can be highly parallelized on clusters
##
## Always check the config.txt file INSIDE the mirrored Singularity directory
## tree, not the one that was initially copied by code/create_dir_tree.sh,
## to set the parameters for this step.
##
## See PHG wiki: https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/Home
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="gvcf-haplos" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=40 #number of cores/tasks
#SBATCH --time=36:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

## Paths to the PHG Singularity image and the base directory created by
## code/create_dir_tree.sh
phg_simg="/project/genolabswheatphg/phg_latest.simg"
base_dir="/project/genolabswheatphg/SRW_test_phg/phg"

## Paths to input data folder and the keyfile
gvcf_in_dir="/project/genolabswheatphg/gvcfs/SRW_single_samp_bw2_excap_GBS_mq20"
keyfile="/project/genolabswheatphg/SRW_test_phg/gVCF_keyfile.txt"


#### Executable ####

echo
echo "Start time:"
date

#### Load Reference Ranges
mkdir -p "$base_dir"
singularity run \
    -B "$base_dir":/tempFileDir/ \
    "$phg_simg" /LoadGenomeIntervals.sh config.txt reference.fa intervals.bed ref_load_data.txt true 


#### Create Haplotypes from gVCFs

## Copy gVCFs and keyfile into mirrored Singularity directory
cp "$keyfile" "$base_dir"/data/keyfile.txt
mkdir "$base_dir"/data/gvcfs
cp "$gvcf_in_dir"/*.vcf.gz "$base_dir"/data/gvcfs/

## Run create haplotypes
singularity run \
    -B "$base_dir":/tempFileDir/ \
    "$phg_simg" /CreateHaplotypesFromGVCF.groovy -config /tempFileDir/data/config.txt

echo
echo "End time:"
date
