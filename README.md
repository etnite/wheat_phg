# wheat_prac_hap_graph

Author: [Brian Ward](https://brianpward.net/)  
License: [GPLv3](https://opensource.org/licenses/GPL-3.0)

## Description

This IS NOT a source code repository for the Practical Haplotype Graph (PHG) 
project. For that, refer to:

    https://bitbucket.org/bucklerlab/practicalhaplotypegraph/src

Rather, this repository is a collection of scripts used to create a practical 
haplotype graph for wheat.

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
can import Docker images, e.g.: 

```
singularity pull docker://maizegenetics/phg
```

However, one issue is that singularity does not appear to
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
