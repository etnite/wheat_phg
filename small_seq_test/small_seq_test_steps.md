# ﻿Small Seq Test

The Small Seq Test simulates data that can take us through the steps of initializing the PHG with a reference sequence and reference ranges, then create haplotypes starting from assembled genomes in fasta, and from non-assembled short reads in fastq.  From there we create consensus haplotypes in the ranges and then submit skim sequence fastq to the PHG to impute variants from skim sequence.

Scripts for the Small Seq Test are in /home/jj332_phg/SmallSeqTest/Scripts
You first want to get them into your own /workdir/<user> folder:
> cd /home/jj332_phg/SmallSeqTest/Scripts
> cp setupDirectories.sh /workdir/<user>  (where <user> is your loginID)
> cd /workdir/<user>
> vim setupDirectories.sh

This is to edit the script to replace “USER=jj332” with your loginID
In vim, hit “i”.  That puts you in “insert” mode.  You can move your cursor with arrows.  Backspace delete “jj332”.  Type in your loginID.  Now hit <esc> to get out of insert mode.  Now type :w to write the changes and :q to quit from vim.

BETTER WAY:

sed -i 's/USER=jj332/USER=bpward2/' *.sh

NOTE: Do NOT change jj332 in the line that begins “RAW_FILES=”

Now run the setupDirectories.sh script:
 ./setupDirectories.sh

That should create two directories in your /workdir/<user>: SSTscripts and SSTtemp.  The PHG will go into SSTtemp.  You should go into SSTscripts to run the scripts from there.

> cd SSTscripts

The first script to run simulates the datasets. createSmallSeqDS.sh script uses a docker to create a small data set for testing.  This data set is created in the /root directory of the docker.  These files are mounted to the specified place in the user's /workdir.  The scripts have docker mount points setup that reflect where the data lives for each step. 

1.  Create the data set:

> ./createSmallSeqDS.sh > createSmallSeqDS.log

Now initialize the PHG with the reference sequence that was simulated and the anchor reference ranges that were chosen.

2.  Load the reference genome data:

> ./loadGenomes.sh > loadGenomes.log

Two “assembled” sequences were simulated.  You can load Sthem and create haplotypes from them using

3.  Load assembly haplotypes:

> ./createAssemblyHaplotypes.sh > createAssemblyHaplotypes.log

Short reads from deep sequencing were also simulated.  You can load them and create haplotypes from them using

4.  Create/load the GATK raw haplotypes:

> ./createHaplotypesAll.sh  > createHaplotypesAll.log

Note that this step takes longer because the PHG software is having to do a lot more aligning than it did to load the assembled sequences.

OK!  Now you have a PHG with a reference sequence and defined reference ranges that are populated with haplotypes, one haplotype per loaded individual.  Within some reference ranges, some individuals may share haplotypes, or nearly so.  To recognize this similarity, you want to collapse haplotypes if they are similar enough.

5.  Create/load consensus data:
  
> ./createConsensus.sh > createConsensus.log

Finally now imagine that you have skim-sequenced a bunch of lines, and you want to use that sequence to infer (impute) the variants carried by the lines.  You submit the skim sequence fastq to the PHG so that it can find the most likely path of the gamete through the PHG.  Having done that, it can output all the variants along that path.

6.  Find and export paths through the graph:

> ./findAndExportPath.sh  > findAndExportPath.log

The variants that are imputed will be exported as a VCF file


In summary, this short pipeline creates and populates an sqlite database.  The data files needed for loading the database are created in the first step (createSmallSeqDS.sh). Starting with loadGenomes.sh, each step above will populate additional data in the sqlite database.  This database is named "phgSmallSeq.db".  The path to the database is:
/workdir/<user>/SSTtemp/phgSmallSeq/phgSmallSeq.db

If you are conversant in SQL, you can query that database.

The simulation step can be tuned to create different datasets that differ in various ways that are interesting for the PHG. For this example, the gene length, inter-gene length, and number of genes can be modified.  In addition, the "refInterGeneDelete" parameter determines structural diversity in the genomes.  The values given to these parameters are set in configSmallSeq.txt  createSmallSeqDS.sh mounts configSmallSeq.txt to a docker directory. Once you have successfully run the pipeline with the default configuration file you may play with the other available parameters.  These exist in the config file but are commented out.

If you wish to re-run createSmallSeqDS.sh a second time, either delete or rename your phgSmallSeq folder.  Docker, however, is unable to delete/overwrite this as it lives in docker in /root.  Thus, you will see an error if you re-run createSmallSeqDS.sh without removing/renaming the old phgSmallSeq directory.  Before deleting files/directories create on cbsu machines, you must reclaim them via the command:
> docker1 claim


Without going into detail, this enables you now to delete the old phgSmallSeq directory.








