#!/usr/bin/env python3

'''
Multinomial Logistic Regression Recursive Feature Elimination

Brian Ward
brian@brianpward.net
https://github.com/etnite

This script is an attempt to create a method to identify a minimal set of 
predictors (SNPs) which fully describe a non-ordinal categorical response 
(haplotypes) for a given region of the genome. My idea is to perform a (possibly 
multinomial) logistic regression, using recursive feature elimination (RFE) to 
select a subset of SNPs with high haplotype prediction accuracy.

There are many different methods for finding "tag SNPs" which can capture the 
variance explained by multiple surrounding SNPs. This problem is non-trivial and 
potentially NP-hard. See https://doi.org/10.1109/CSB.2005.22 for one example and 
some background on different methods. Much of this work was performed 
around the time of the first human HapMap project (ca. 2005), and hence there are 
many links to software on websites which no longer exist. These programs are 
also likely to be written in C or Perl, and getting them to run and use current 
data formats can be very difficult.

Many tag SNP identification algorithms are designed to find sets of tag SNPs 
across chromosomes/genomes, without any phenotypic data. I think that my particular
use-case here is actually somewhat simpler, since I:

1. Have access to a "phenotype" (haplotype groups determined by a simple distance 
   matrix using all SNPs/variants in the region of interest)
2. Only need to focus on a single region, obviating the need to find haploblock 
   boundaries

This script performs repeated k-fold validation for a user-specified range of
SNP numbers. For instance, one could run 10 repeats of 5-fold validation using
from 1 to 10 SNPs. The script will also always run the k-fold validation using
all SNPs within the specified region.

scikit-learn's RFE algorithm is used for feature elimination, and we use the
handy scikit-allel module to import data from a VCF file: 
http://alimanfoo.github.io/2017/06/14/read-vcf.html
'''

import os
import sys
import random
import allel
import numpy as np
import pandas as pd
from sklearn.feature_selection import RFE
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, KFold
from sklearn import preprocessing
import seaborn as sns
random.seed(123)


#### User-Defined constants ####

## Path to input VCF file, and region specified as CHROM:START-END
vcf_file = '/home/brian/repos/manuscripts/manu_2018_stripe_rust/input_data/geno/GAWN_KM_Yr_postimp_filt.vcf.gz'
reg = '4B:5000000-7000000'

## Haplotypes file (.csv with header and at least two columns - 'sample' and 
## the specified response column)
haps_file = '/home/brian/Downloads/haplotype_groupings.csv'
response = '3BS'

## Path to directory to save output files
out_dir = '/home/brian/Downloads/3BS_mlogit_null_test_out'

## Cross-val parameters - range of SNP numbers, repeats, and number of folds
min_snps = 1
max_snps = 10
n_reps = 10
val_k = 5


#### Executable ####

os.makedirs(out_dir, exist_ok = True)

## Read in the VCF file using scikit-allel
vcf = allel.read_vcf(vcf_file, region = reg)
if vcf is None:
    sys.exit("\nNo variants present in specified region!")
print("\nNumber of Samples:", len(vcf['samples']))
gt = allel.GenotypeArray(vcf['calldata/GT'])

## Adjust max_snps if necessary
preds = vcf['variants/ID']
print("Number of SNPs in region:", str(len(preds)))
if max_snps > len(preds):
    max_snps = len(preds)
    print("\nWARNING - max_snps was set to a value higher than the number of SNPs in the region\n\nAdjusting max_snps to", str(len(preds)))


## Convert the genotypic matrix to dataframe in minor allele dosage format
dos = gt.to_n_alt()
dos = dos.transpose()
dos_df = pd.DataFrame(dos)
dos_df.columns = vcf["variants/ID"]
dos_df["sample"] = vcf["samples"]


## Read in the line haplotype information
haps = pd.read_csv(haps_file)
haps = haps[["sample", response]]


## Merge together predictors (SNPs) and response (haplotype groupings)
## Set response vector and predictors matrix
merged = pd.merge(haps, dos_df, how = "inner", on = "sample")
print("\nNumber of samples matching between VCF and haplotypes file:", str(merged.shape[0]))
y = merged[response]
X = merged[vcf["variants/ID"]]


## Standardize SNPs
sc = preprocessing.StandardScaler()
X_std = sc.fit_transform(X)


## Initialize Dataframe for output
d = {'rep': list(range(1, n_reps + 1)) * val_k,
     'fold': list(range(1, val_k + 1)) * n_reps}
d['rep'].sort()
proto_df = pd.DataFrame(data = d)


## Create multinomial logistic regressor
mlogit = LogisticRegression(solver = 'lbfgs', penalty = 'none', 
    multi_class = 'auto', max_iter = 5e3)


## Loop through number of SNPs
max_snps = max_snps + 1
cv_outlist = [0] * len(range(min_snps, max_snps))
for i, n_snps in enumerate(range(min_snps, max_snps)):
    sub_score_arr = np.zeros((n_reps, val_k), dtype = np.float)

    ## Perform recursive feature elimination; subset SNP matrix
    rfe = RFE(mlogit, n_snps)
    fit = rfe.fit(X_std, y)
    pred_str = ' + '.join(preds[fit.support_])
    X_sub = X_std[:, fit.support_]

    ## Run repeated cross-validation
    for j in range(n_reps):
        k_fold = KFold(n_splits = val_k, random_state = j, shuffle = True)
        sub_score_arr[j,:] = cross_val_score(mlogit, X_sub, y, cv = k_fold, scoring = 'accuracy')

    ## Record cross-validation results
    cv_outlist[i] = proto_df.copy()
    cv_outlist[i]['nSNPs'] = n_snps
    cv_outlist[i]['accuracy'] = sub_score_arr.flatten()
    cv_outlist[i]['SNP_IDs'] = pred_str

## Concatenate the DFs together
cv_concat = pd.concat(cv_outlist)


## One final iteration using all SNPs in region
## iff max_snps < number of all SNPs
if max_snps < len(preds):
    sub_score_arr = np.zeros((n_reps, val_k), dtype = np.float)
    for j in range(n_reps):
        k_fold = KFold(n_splits = val_k, random_state = j, shuffle = True)
        sub_score_arr[j,:] = cross_val_score(mlogit, X_std, y, cv = k_fold, scoring = 'accuracy')
    proto_df['nSNPs'] = X_std.shape[1]
    proto_df['accuracy'] = sub_score_arr.flatten()
    proto_df['SNP_IDs'] = 'all'
    cv_concat = pd.concat([cv_concat, proto_df])


## Calculate summary stats
cv_concat = cv_concat[['nSNPs', 'rep', 'fold', 'accuracy', 'SNP_IDs']]
cv_summ = cv_concat[['nSNPs', 'accuracy']].groupby('nSNPs').describe()
cv_summ.columns = cv_summ.columns.droplevel(0)
cv_sem = cv_concat[['nSNPs', 'accuracy']].groupby('nSNPs').sem()
cv_sem = cv_sem.rename(columns = {'accuracy': 'se'})
cv_summ = pd.merge(cv_summ, cv_sem, how = "inner", on = 'nSNPs')
cv_summ = cv_summ[['count', 'mean', 'std', 'se']]


## Create barplot
sns.set()
acc_plot = sns.barplot(x = 'nSNPs', y = 'accuracy', data = cv_concat, color = 'royalblue')
figure = acc_plot.get_figure()
figure.set_size_inches(10, 8)


## Output files and figures
sub_preds = list(preds[fit.support_])
merged[['sample', response] + sub_preds].to_csv(os.path.join(out_dir, 'haplotypes_SNPs.csv'), na_rep = 'NA', index = False)
cv_concat.to_csv(os.path.join(out_dir, 'CV_accuracies.csv'), na_rep = 'NA', index = False)
cv_summ.to_csv(os.path.join(out_dir, 'CV_summaries.csv'), na_rep = 'NA', index = True)
figure.savefig(os.path.join(out_dir, 'CV_barplot.png'), dpi = 72)
