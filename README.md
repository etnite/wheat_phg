# wheat_prac_hap_graph

Author: [Brian Ward](https://brianpward.net/)  
License: [GPLv3](https://opensource.org/licenses/GPL-3.0)

## Description

This IS NOT a source code repository for the Practical Haplotype Graph (PHG) 
project. For that, refer to:

    https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src

Rather, this repository is a collection of scripts used to create a practical 
haplotype graph for wheat.

**NOTE:** All scripts ending in "parallel.sh" are intended to be run using the
code/arrayer.sh helper script, which is used to submit arrays of independent jobs
on a cluster. 

Code in this repository may be used for either:

1. Building a practical haplotype graph
2. Performing the steps necessary for classical short variant calling
using Illumina sequencing data.

For either of these two tasks, a typical initial workflow would be:

1. run code/concat_fastqs.sh to generate a single mated pair of fastq files for
each individual sample
2. run code/fastq_filt_trim/bbduk_filt_trim_paired_parallel.sh to clean fastq reads
3. run code/alignment/bowtie2_align_parallel.sh to align reads using bowtie2 and create
single sample bam alignment files
4. run code/filt_bam_files_parallel.sh to filter bam files after alignment

For subsequent loading of data into a PHG, the user can then optionally run
code/create_gvcf.sh

For traditional variant calling, the following steps may be performed after bam file
generation:

1. run code/call_variants/call_variants.sh
2. run code/call_variants/filter_raw_vcf.sh
3. run code/call_variants/rename_annotate_split_vcf.sh to (optionally) rename samples,
predict variant transcriptional effects using a genome annotation file, and create subset 
VCF files consisting of just indels and just SNPs.

These steps should create VCF files suitable for subsequent imputation (if necessary).

## Dependencies

### Operating system and hardware

Many scripts in this repository assume that the user is using a high-performance
computing cluster running a 64-bit Linux distribution. It is assumed that the cluster uses
SLURM scheduling.

### Docker or Singularity

Although the PHG can be installed from source, it is designed to be used through
a docker image to handle all software dependencies and achieve platform
independence: 

```
docker pull maizegenetics/phg
docker pull maizegenetics/phg_postgres
```

Note that many clusters use Singularity instead of Docker. Singularity
can import Docker images. However, one issue is that singularity does not appear to
support the creation of nested directory structures, which the PHG requires. Therefore,
**ALL DIRECTORIES USED WITH SINGULARITY MUST BE MANUALLY SPECIFIED**

More to come on this later...

## small_seq_test

This directory holds a set of scripts and data for running a small simulated PHG
creation. Things to note:

* The lines in all the scripts that start with "PHG_DIR=" will need to be modified
in order to point to the correct directory. Currently the method for doing this is
a bit circuitous, so I may modify in the future.
* The scripts currently use the command docker1, which is specific to the Cornell
BioHPC. If running on another computer/cluster, these commands should be changed to
just "docker"
