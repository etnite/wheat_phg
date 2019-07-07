#!/bin/bash

## Paired end alignment using Bowtie2
##
## This script aligns a single pair of mated fastq files.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for corresponding fastq files using the pattern:
## fastq_dir/samplename_R[12].fastq[.gz]
##
## The output consists of a single sorted BAM file, with PCR duplicates marked.
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
##   5) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined constants ####

## Reference genome fasta ("ref") must already be indexed using bowtie2-build
## and samtools index
fastq_dir="/project/genolabswheatphg/test/BPW_pipeline_test/filt_fastq"
out_dir="/project/genolabswheatphg/test/BPW_pipeline_test/alignments"
ref="/project/genolabswheatphg/v1_refseq/Clay_splitchroms_reference/161010_Chinese_Spring_v1.0_pseudomolecules_parts.fasta"


#### Executable  ####

module load bowtie2
module load samtools

echo
echo "Start bowtie2_align_parallel.sh"
echo "Start time:"
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
fq1="${fastq_dir}"/"${samp}"_10K_R1.fastq.gz
fq2="${fastq_dir}"/"${samp}"_10K_R2.fastq.gz

## For some reason, the barcode indexes in FASTQ files can contain some N
## values for the first few reads. I don't know the significance of this. 
## Let's just grab line 10,001:
## NOTE: In the case of concatenated fastq files, flowcells, lanes,
## barcodes, etc. may all vary. Therefore these sections are disabled
#one_line=$(zcat $fq1 | head -n 10001 | tail -n -1)

## This sets the first four fields (instrument:run_id:flowcell:lane) of the 
## first line of the fastq file as a variable named "id"
#id=$(echo "$one_line" | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')

## Now get the barcode from the 10th field of the line
## NO
#bar=$(echo "$one_line" | cut -f 10 -d":" | sed 's/+/-/')


ref="${ref%.*}"
bowtie2 -x "${ref}" \
        --threads $SLURM_NTASKS \
        --rg-id "${samp}" \
        --rg SM:"${samp}" \
        --rg PL:ILLUMINA \
        --sensitive-local \
        --phred33 \
        -1 "${fq1}" \
        -2 "${fq2}" |
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

echo
echo "End time:"
date
