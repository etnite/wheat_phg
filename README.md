# wheat_prac_hap_graph

Author: [Brian Ward](https://brianpward.net/)  
License: [GPLv3](https://opensource.org/licenses/GPL-3.0)

## Description

This IS NOT a source code repository for the Practical Haplotype Graph (PHG) 
project. For that, refer to:

[https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src](https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src)

Rather, this repository is a collection of scripts used to create a practical 
haplotype graph for wheat.

## Dependencies

### Operating system and hardware

Many scripts in this repository assume that the user is using a high-performance
computing cluster running a 64-bit Linux distribution. It is assumed that the cluster uses
SLURM scheduling.

### Docker and Singularity

Although the PHG can be installed from source, it is designed to be used through
a docker image to handle all software dependencies and achieve platform
independence: 

```
docker pull maizegenetics/phg
docker pull maizegenetics/phg_postgres
```

Note that many clusters use Singularity instead of Docker. Singularity
can import Docker images, e.g. to create a Singularity image called "wheat.simg",
run: 

```
singularity build  wheat.simg docker://maizegenetics/phg:0.0.17
```

To then look inside the created Singularity image, you can run:

```
singularity shell wheat.simg
cd /
```

`singularity shell` will run a subshell within the Singularity image.
However, the working directory remains the same, so the "cd /" is required to
access the image's root directory.

The main issue with using Singularity instead of Docker is that Singularity does
not run with the privaleges necessary to create new directories. Therefore,
steps in which the PHG Docker scripts create directories will fail when using
Singularity. **code/create_dir_tree.sh** is intended to create an external directory
structure which mirrors the internal directory structure used by the Docker scripts,
thereby circumventing the directory creation problems.


## small_seq_test

This directory holds a set of scripts and data for running a small simulated PHG
creation. Things to note:

* The lines in all the scripts that start with "PHG_DIR=" will need to be modified
in order to point to the correct directory. Currently the method for doing this is
a bit circuitous, so I may modify in the future.
* The scripts currently use the command docker1, which is specific to the Cornell
BioHPC. If running on another computer/cluster, these commands should be changed to
just "docker"
