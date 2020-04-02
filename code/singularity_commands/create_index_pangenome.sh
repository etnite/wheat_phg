#! /bin/bash

## Create PHG Pangenome
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This wrapper script will proceed through the steps necessary to create and
## index a PHG pangenome fasta file. If the user supplies a consensus method
## name (parameter consens_meth), then Consensi will be created prior to indexing
## a pangenome file. However, see description in user-defined constants section
## below for more notes on different possibilities for creating haplotypes and/or
## consensi.
##
## Always check the config.txt file INSIDE the mirrored Singularity directory
## tree, not the one that was initially copied by code/create_dir_tree.sh,
## to set the parameters for this step.
##
## See PHG wiki: https://bitbucket.org/bucklerlab/practicalhaplotypegraph/wiki/Home
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="pangen" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=40 #number of cores/tasks
#SBATCH --time=4:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

## Paths to the PHG Singularity image and the base directory created by
## code/create_dir_tree.sh
phg_simg="/project/genolabswheatphg/phg_latest.sif"
base_dir="/project/genolabswheatphg/SRW_1ABD_phg_test/phg"

## Paths to input data folder and the keyfile
gvcf_in_dir="/project/genolabswheatphg/gvcfs/SRW_single_samp_bw2_excap_GBS_mq20"
keyfile="/project/genolabswheatphg/SRW_1ABD_phg_test/gVCF_keyfile.txt"

## Set names of haplotype and consensus method
## Haplotype method is generally set to "GATK_PIPELINE"
## If consens_meth == haplo_meth, then consensus building is skipped,
##   and pangenome is built with all haplotypes
##   Otherwise consens_meth can be any convenient name for the consensus method used
## Can get more creative - for instance load in reference sequence as a haplotype
## using <haplotype_method>:<reference_method>
## use select * from methods; on PHG SQLite DB to see the reference method's name
haplo_meth="GATK_PIPELINE"
consens_meth="CONSENSUS_MXDIV_1E-7"

## The number of bases for minimap2 to load into DB for each batch in format, e.g. "100G"
## Typically set to > number of bases in pangenome
## Rough estimate of pangenome size is (portion genome in ref. ranges x number taxa)
## Though number will be much lower if consensus is run
n_base_load="30G"


#### Executable ####

echo
echo "Start time:"
date

## create consensus indicator variable
run_consens=$(echo "$consens_meth" | tr '[:upper:]' '[:lower:]')

## Pangenome creation/indexing
if [[ "$consens_meth" != "$haplo_meth" ]]; then

    ## Create Consensus
    singularity run \
        -B "$base_dir":/tempFileDir/ \
        "$phg_simg" /CreateConsensi.sh /tempFileDir/data/config.txt reference.fa "$haplo_meth" "$consens_meth" 

fi

## Index the pangenome
    singularity run \
        -B "$base_dir":/tempFileDir/ \
        "$phg_simg" /IndexPangenome.sh pangenome config.txt "$consens_meth" "$n_base_load" 21 11x

echo
echo "End time:"
date
