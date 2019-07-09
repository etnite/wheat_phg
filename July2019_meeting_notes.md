# PHG July 2019 Meeting

## Monday, 08July2019

### Accessing BioHPC - ssh and FTP

my rental machine: cbsumm30
Jason's rental machine: cbsumm31

Filezilla FTP connection (port 22):
sftp://cbsumm30.biohpc.cornell.edu
sftp://cbsulogin.biohpc.cornell.edu

Corresponding ssh connection:
ssh bpward2@cbsumm30.biohpc.cornell.edu

Home directory:

	/home/bpward2

Work directory:

	/workdir/bpward2

Outside of Cornell (after workshop for a while), can access the rental machines
by logging in to the general login nodes:

* cbsulogin

And then ssh to the rental machines


### What is the Practical Haplotype Graph?

Nomenclature:

* reference ranges are... ranges along the reference genome, some of which will
be used as *anchors* (see below), which are *generally* conserved/easy to align
(but which still contain some variation)
* reference range = reference interval = genome interval
* anchor = genic interval

Note that database may contain input haplotypes and consensus halpotypes.
Creation of consensus haplotypes is tunable and optional.

In maize, there is about 10x variation between haplotypes, vs. within haplotypes.
This ratio varies by species.

Biparental mapping: can be restricted to two parent lines to find most likely
path through haplotypes for progeny. However, it cannot yet make use of 
*a priori* estimates of recombination (i.e. marker genetic positions)

* "Gamete" in the PHG is basically synonymous with "genotype" or "line", as we
are dealing only with inbred lines.

* For assemblies, Buckler lab uses Mummer4 for whole-genome alignment (e.g.
identification of inversions, translocations, insertions, deletions, etc.)