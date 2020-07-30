## Create a keyfile for gVCF files


#### User-Defined Input/Output #################################################

## The gVCF list file would typically be created with the Linux find or realpath 
## command, e.g. in the folder containing all the gVCF files run:
## realpath *.g.vcf.gz | sort > gVCF_list.txt
gvcfs_list_file <- "gvcf_list.txt"

## Path to output keyfile - this will be a tab-delimited text file
out_key_file <- "gvcf_keyfile.txt"


#### Create the keyfile - some hand modification may be necessary below ########

gvcfs_list <- scan(gvcfs_list_file, what = "character", sep = "\n")

key <- data.frame("files" = gvcfs_list)

## Get gVCF sample names
key$sample_name <- basename(gvcfs_list)
key$sample_name <- sub(".g.vcf.gz", "", key$sample_name)

## Here we're just setting sample_description = sample_name
key$sample_description <- key$sample_name

## Set some constant columns
key$type <- "GVCF"
key$chrPhased <- "true"
key$genePhased <- "true"

## NOTE - In general this value should be close to 1 for inbred crops. However,
## I don't fully know how to find an ideal value
key$phasingConf <- 0.95

## Sort columns and output
key <- key[c("sample_name", "sample_description", "files", "type", "chrPhased", 
             "genePhased", "phasingConf")]
write.table(key, file = out_key_file, sep = "\t", row.names = FALSE, quote = FALSE)