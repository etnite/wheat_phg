#!/bin/bash

## Paired end alignment using Bowtie2
##
## This script aligns a single pair of mated fastq files.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for corresponding fastq files using the pattern:
## fastq_dir/samplename_R[12].fastq[.gz]
##
## The output consists of a single sorted BAM file.
##
## This script is intended to be used with parallelize.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##
##   1) bowtie2 is run using the --sensitive-local option. This can be changed
##      manually in the call to bowtie2 below
##   2) The reference genome fasta must be indexed using samtools faidx
##   3) The script will produce either a .bai or a .csi index of the output
##      BAM file, depending on the length of the largest chromosome/contig in
##      the reference fasta
##   4) This script requires two sorting steps in samtools, so can be a bit
##      slow...
################################################################################


#### User-defined constants ####

## Reference genome fasta ("ref") must already be indexed using bowtie2-build
fastq_dir="/project/genolabswheatphg/raw_data/wheatCAP_parents"
out_dir="/project/genolabswheatphg/filt_fastqs/wheatCAP_parents"
ref=""


#### Executable  ####

module load bowtie2
module load samtools

date
mkdir -p "${out_dir}"
samp=$1

## Check if reference genome fasta index exists
if [[ ! -f "${ref}".fai ]]; then
    echo "Reference genome must be indexed using samtools faidx command"
    exit 1;
fi

## Get size of largest chromosome/contig in reference
max_chr=$(cut -f2 "${ref}".fai | sort -nr | head -n 1)

## Set forward and reverse read fastq files
fq1="${fastq_dir}"/"${samp}"_R1.fastq*
fq2="${fastq_dir}"/"${samp}"_R2.fastq*

## For some reason, the barcode indexes in FASTQ files can contain some N
## values for the first few reads. I don't know the significance of this. 
## Let's just grab line 10,001:
one_line=$(zcat $fq1 | head -n 10001 | tail -n -1)

## This sets the first four fields (instrument:run_id:flowcell:lane) of the 
## first line of the fastq file as a variable named "id"
id=$(echo "$one_line" | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')

## Now get the barcode from the 10th field of the line
bar=$(echo "$one_line" | cut -f 10 -d":" | sed 's/+/-/')

bowtie2 -x "${ref}" \
        -1 "${fq1}" \
        -2 "${fq2}" \
        --rg-id $(echo "@RG\tID:${id}\tSM:${i}\tLB:${id}_${i}\tBC:${bar}\tPL:ILLUMINA") \
        --phred33 \ 
        --sensitive-local |
        samtools sort -n -O SAM - |
        samtools fixmate -m -O SAM - - |
        samtools sort -O SAM - |
        samtools markdup - "${out_dir}"/"${samp}".bam

## Default .bai indices can only handle contigs up to 2^29 bases. If any contig
## exceeds this length, use the more robust .csi index
if [[ $max_chr > 536870912 ]]; then
    samtools index -c "${out_dir}"/"${samp}".bam
else
    samtools index -b "${out_dir}"/"${samp}".bam
fi

date
