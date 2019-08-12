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

## Wednesday, 10July2019

Filtering gVCFs - can filter on QUAL (overall base quality??), and on GQ 
(genotype quality)


Send AGS2000 data exome capture to:

dllarkin@uark.edu

### Reference Ranges - Dan Ilut

Good reference range features:

1. Consistently genotyped across samples (e.g. exome capture sites)
2. Contain variability among lines
3. Distinct among the wheat sub-genomes (confounded with #2)

Partly depends upon intended use - e.g. genomic selection vs. gene cloning. For 
genomic selection, it doesn't matter as much if anchors are "tuned" exceptionally
well

For gene cloning, must take into account the _region of interest_. There
shouldn't be huge jumps in the number of haplotypes between adjacent anchors (though
depends on distance between anchors).

Missing data shouldn't matter as long as there is enough data present to define
different haplotypes.

May be good idea to filter anchors by missing data levels.

First build PHG without any collapsing (or with very strict collapsing thresholds
to just remove sequencing error), evaluate, and then implement collapsing
as desired.

For building PHG, use the cleanest, highest-confidence data possible (err on the
side of quality over quantity)

Need to quantify divergence between subgenomes (within variety) vs. divergence 
between varieties (at the same locus)


## Thursday, 11JULY2019

### Conference Call


