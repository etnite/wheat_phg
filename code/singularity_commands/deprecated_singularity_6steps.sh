#! /bin/bash

#### SLURM job control #### 

#SBATCH --job-name="phg-singularity" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=40 #number of cores/tasks
#SBATCH --time=10:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

phg_simg="/project/genolabswheatphg/phg_latest.sif"
base_dir="/project/genolabswheatphg/SRW_1ABD_phg_test/phg"
gvcf_in_dir="/project/genolabswheatphg/gvcfs/SRW_single_samp_bw2_excap_GBS_mq20"
gvcf_keyfile="/project/genolabswheatphg/SRW_1ABD_phg_test/gVCF_keyfile.txt"

haplo_meth="GATK_PIPELINE"
consens_meth="CONSENSUS_MXDIV0001"

#### Executable ####

echo
echo "Start time:"
date

#### Load Reference Ranges
#singularity run \
#    -B "$base_dir":/tempFileDir/ \
#    "$phg_simg" /LoadGenomeIntervals.sh config.txt reference.fa intervals.bed ref_load_data.txt true 


#### Create Haplotypes from gVCFs

## Copy gVCFs and keyfile into mirrored Singularity directory
#cp "$gvcf_keyfile" "$base_dir"/data/keyfile.txt
#cp "$gvcf_in_dir"/*.vcf.gz "$base_dir"/data/gvcfs/

## Run create haplotypes
#singularity run \
#    -B "$base_dir":/tempFileDir/ \
#    "$phg_simg" /CreateHaplotypesFromGVCF.groovy -config /tempFileDir/data/config.txt


#### Create Consensi

#singularity run \
#    -B "$base_dir":/tempFileDir/ \
#    "$phg_simg" /CreateConsensi.sh /tempFileDir/data/config.txt reference.fa "$haplo_meth" "$consens_meth" 

#### Index PanGenome
singularity run \
    -B "$base_dir":/tempFileDir/ \
    "$phg_simg" /IndexPangenome.sh pangenome config.txt "$consens_meth" 30G 21 11x 

#singularity run \
#-B /home/jason.fiedler/phg/HRSW_test/tempFileDir/:/tempFileDir/ \
#/home/jason.fiedler/phg/phg_latest.sif \
#/FindPathMinimap2.sh pangenome500 config.txt \
#GATK_PIPELINE GATK_PIPELINE \
#Dummy PATH3 \
#/tempFileDir/data/fastq/genotypingKeyFile500SR.txt \
#/tempFileDir/data/fastq/genotypingKeyFile_pathKeyFile500SR.txt 
#
#singularity run \
#-B /home/jason.fiedler/phg/HRSW_test/tempFileDir/:/tempFileDir/ \
#/home/jason.fiedler/phg/phg_latest.sif \
#/ExportPath.sh config.txt GATK_PIPELINE 500.path3out PATH3 

echo
echo "End time:"
date
exit 0;
