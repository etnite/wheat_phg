#!/bin/bash -x
# This script sets up directories for the Small Seq Test (SST)
# Edit the script to replace jj332 with your user id
# The variable USER is usually set to your user id, but just in case
USER=bpward2

RAW_FILES=/home/jj332_phg/SmallSeqTest/Scripts

# Move files from the phg raw files directory to your working directory
SCRIPT_DIR=/workdir/${USER}/SSTscripts/
mkdir -p $SCRIPT_DIR
chmod 777 $SCRIPT_DIR
cp ${RAW_FILES}/* ${SCRIPT_DIR}

# Create a directory for PHG work
PHG_DIR=/workdir/${USER}/SSTtemp/
mkdir -p $PHG_DIR
chmod 777 $PHG_DIR
cp ${RAW_FILES}/configSmallSeq.txt ${PHG_DIR}
