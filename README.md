# wheat_prac_hap_graph

## Description

This IS NOT a source code repository for the Practical Haplotype Graph (PHG) 
project. For that, refer to:

    https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src

Rather, this repository is a collection of scripts used to create a practical 
haplotype graph for wheat. The first step is to set up directories and input
data. This repository does not house the required input data, as much of
this data consists of large files. Instead, it sets up a separate directory
structure where data is downloaded and then processed.

## Dependencies

### Operating system and hardware

For the actual graph construction, a 64-bit Linux distribution is required,
along with beefy computational resources (e.g. multiple cores and a significant
amount of RAM). Subsequent use of the graph to call/impute genotypes should
be much less resource-intensive.

### Docker

Although the PHG can be installed from source, it is designed to be used through
a docker image to handle all software dependencies and achieve platform
independence. If you are working on a server or HPC cluster, Docker may already 
be installed. Otherwise Docker can be installed and maintained with the package 
management systems of several popular Linux distributions:

* Ubuntu
* Debian
* CentOS
* Fedora

and can be installed from source on other distributions. Instructions below 
are for installing Docker Community Edition (CE) on an Ubuntu 16.04 distribution.
See https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce
for more details:

```
## Update apt cache
$ sudo apt-get update

## Install packages to allow apt to use repositories over HTTPS
$ sudo apt-get install \
    apt-transport-https \
    ca-certificates \
    curl \
    software-properties-common

## Download Docker's official GPG key
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -

## check that key has fingerprint 9DC8 5822 9FC7 DD38 854A E2D8 8D81 803C 0EBF CD88
$ sudo apt-key fingerprint 0EBFCD88

## Add the stable repository
$ sudo add-apt-repository \
   "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
   $(lsb_release -cs) \
   stable"

## Update apt cache and install docker
$ sudo apt-get update
$ sudo apt-get install docker-ce
```

Once docker is installed, the PHG image can be installed by running the following
command in a terminal:

```
docker pull maizegenetics/phg
docker pull maizegenetics/phg_postgres
```

## Directory setup and download

The script dir_data_setup.sh will attempt to set up a directory and data
structure that is suitable for running the PHG. Note that this script will
require some modification over time as data sources change.

## Config directory

The directory config/ holds some sample configuration files from the Buckler
lab. When dir_data_setup.sh is run it will copy these files into the new
project directory - that is, to [project_directory]/config. Note that these files
should be modified as necessary after copying to customize the graph construction.


## Steps in PHG Creation

The creation of a PHG has four main steps:

1. Identify reference intervals (usually genic regions) for inclusion in the
   graph. 
2. Load the reference genome
3. Create haplotypes
4. Create consensus haplotypes


## Identifying reference intervals

This step uses an input reference genome fasta file and its corresponding
annotation (i.e. .gff3) file to identify genes. This step only identifies the
position of genic regions, so overlapping genes will be collapsed into a single
region.

This step can optionally expand reference intervals outwards from genes until
repetitive sequence is encountered. **This is not necessary if working with
exome capture data**. To do this, a kmer-analysis is performed to search for the
the top repetitive kmers in the genome. The critical parameters for performing
this analysis are *e* - the number of bases by which to expand the genic regions
for initial interval selection, *p* - the proportion of kmers considered
repetitive (e.g. the most frequent 5% of all kmers), *n* - the number of copies
over which a kmer is considered repetitive (overides -p).

The output consists of a .bed file identifying the genic regions, as well as the
expanded genic regions.

Note that for exome capture data, which should theoretically consist of only
genic regions, -e should simply be set to 0. The script 
docker_run_commands/create_ref_intervals.sh is designed to perform this step of
the process with wheat exome capture data.

**NOTE**: One wheat-specific detail - some chroms are larger than 512Mb, which is
a problem for short-read aligners. The IWGSC v1.0 RefSeq has each chromosome
split into two halves. The way they split them was not very elegant -
chromosome 6D appears to have been split inside a gene, so the gene is present in 
both halves of the chrom. The start position of the gene is 450509070, so it can
be grepped out of the .bed file before proceeding (might have to add 1 to the
position - I can never remember). Many of the wheat chroms are less than 512Mb,
so the halves of these chroms could be appended together, but maybe that's more
trouble than it's worth.