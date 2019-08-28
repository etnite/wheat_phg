#!/bin/bash
set -e
#source /home/gbg_lab_admin/miniconda3/bin/activate bwa_align_call


## Filter VCF file - single threaded
##
## NOTE: The max_miss parameter works in an opposite manner of the VCFTools
## implementation, which I find unintuitive. Here, if max_miss is set to 1, SNPs with 100%
## missing data would theoretically be allowed, while setting max_miss to 0 will only
## allow SNPs without any missing data
##
## snpgap will remove SNPs within n bases of an indel (or overlapping)
## indelgap will thin clusters of indels within n bases of each other, to only retain one
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="filt-vcf" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=06:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

vcf_in="/home/gbg_lab_admin/Array_60TB/Wheat_GBS/GAWN_and_KM_Yr/GAWN_KM_50miss_filt/prefilt_and_imp/GAWN_KM_Yr_50miss_prefilt.vcf.gz"
vcf_out="/home/gbg_lab_admin/Downloads/bcftools_filt_test/bcftools_filt_test.vcf.gz"
taxa_list="/mnt/DATAPART2/brian/repos/manuscripts/manu_2019_stripe_rust/entry_lists_and_randomization/KM_entries.txt"
min_maf=0.1
max_miss=0.2
max_het=0.1
min_dp=0
max_dp=1e6
remove_unal="false"
snpgap=3
indelgap=3


#### Executable #####

module load bcftools
#module load vcftools

echo
echo "Start time:"
date

## Grab first letter of remove_unal
remove_unal=${remove_unal:0:1}

out_dir=$(dirname "${vcf_out}")
mkdir -p "${out_dir}"
temp_dir="$(mktemp -d -p "${out_dir}")"
if [[ ! "${temp_dir}" || ! -d "${temp_dir}" ]]; then
    echo "Could not create temporary directory"
    exit 1;
fi

## Echo input parameters to output dir
echo -e "Input VCF\t${vcf_in}" > "${out_dir}"/filtering_params.txt
echo -e "Output VCF\t${vcf_out}" >> "${out_dir}"/filtering_params.txt
echo -e "Taxa subset list\t${taxa_list}" >> "${out_dir}"/filtering_params.txt
echo -e "Minimum MAF\t${min_maf}" >> "${out_dir}"/filtering_params.txt
echo -e "Max missing proportion\t${max_miss}" >> "${out_dir}"/filtering_params.txt
echo -e "Max het. proportion\t${max_het}" >> "${out_dir}"/filtering_params.txt
echo -e "Min. average depth\t${min_dp}" >> "${out_dir}"/filtering_params.txt
echo -e "Max average depth\t${max_dp}" >> "${out_dir}"/filtering_params.txt
echo -e "Unaligned contigs removed?\t${remove_unal}" >> "${out_dir}"/filtering_params.txt
echo -e "SNP overlap gap\t${snpgap}" >> "${out_dir}"/filtering_params.txt
echo -e "Indel overlap gap\t${indelgap}" >> "${out_dir}"/filtering_params.txt

## If taxa_list exists, use to subset samples
## Otherwise retain all samples present in VCF file
if [[ -f $taxa_list ]]; then
	cp $taxa_list "${temp_dir}"/taxa_list.txt
else
	bcftools query --list-samples $vcf_in > "${temp_dir}"/taxa_list.txt
fi

## Some real bcftools power-user stuff here
echo "Filtering VCF..."
echo
if [[ $remove_unal == [Tt] ]]; then
    bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools view - \
        --targets ^UN \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || AVG(FORMAT/DP) < ${min_dp} || AVG(FORMAT/DP) > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het} " \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type z \
        --output "${vcf_out}"
elif [[ $remove_unal == [Ff] ]]; then
	bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools view - \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || AVG(FORMAT/DP) < ${min_dp} || AVG(FORMAT/DP) > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het} " \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type z \
        --output "${vcf_out}"
else
	echo "Please supply 'true' or 'false' for remove_unal"
	exit 1;
fi

bcftools index -c "${vcf_out}"


## Generate summary stats using TASSEL
#echo "Generating summary statistics with TASSEL..."
#$TASSEL_PL -vcf "${vcf_out}" \
#           -GenotypeSummaryPlugin \
#           -endPlugin \
#           -export "${out_dir}"/summary

rm -rf $temp_dir
#source deactivate

echo
echo "End time:"
date
