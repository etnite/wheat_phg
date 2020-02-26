#!/bin/bash

## Very Simple Script to create a Singularity Image from PHG Docker Image
#########################################################################


#### User-defined Constants ####

## Path to output Singularity image. Should have ".simg" extension
simg_path=""


#### Executable ####

singularity build  wheat.simg docker://maizegenetics/phg:0.0.17


## If you want to look inside the image:
## Only run this interactively:
#singularity shell <image_name.simg>
#cd /
