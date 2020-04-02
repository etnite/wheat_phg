#!/usr/bin/env python3

'''
Brian Ward
brian@brianpward.net
https://github.com/etnite
'''

from sklearn.feature_selection import RFE
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, KFold
from sklearn import preprocessing
import allel
import numpy as np
import pandas as pd


## User-Defined constants
vcf_file = '/home/brian/repos/manuscripts/manu_2018_stripe_rust/input_data/geno/GAWN_Yr_postimp_filt_KASP.vcf.gz'
reg = '4B:526943215-598847646'
haps_file = '/home/brian/Downloads/haplotype_groupings.csv'
response = '4BL'
min_snps = 4
max_snps = 5
n_reps = 5
val_k = 5


## Read in the VCF file

vcf = allel.read_vcf(vcf_file, region = reg)
preds = vcf['variants/ID']
gt = allel.GenotypeArray(vcf['calldata/GT'])
#print("Number of Samples:", len(vcf['samples']))
#print("Number of SNPs:", len(preds))
#gt



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
y = merged[response]
X = merged[vcf["variants/ID"]]


## Standardize SNPs

sc = preprocessing.StandardScaler()
X_std = sc.fit_transform(X)
#X_std

## Initialize Output Structures
d = {'rep': list(range(1, n_reps + 1)) * val_k,
     'fold': list(range(1, val_k + 1)) * n_reps}
d['rep'].sort()
proto_df = pd.DataFrame(data = d)



## Multinomial logistic regression
mlogit = LogisticRegression(solver = 'lbfgs', penalty = 'none', multi_class = 'auto', max_iter = 5e3)









#print("Retained features array first 5 rows:")
#X_sub[:5]


max_snps = max_snps + 1
cv_list = [0] * len(range(min_snps, max_snps))
for i, n_snps in enumerate(range(min_snps, max_snps)):
    sub_score_arr = np.zeros((n_reps, val_k), dtype = np.float)

    rfe = RFE(mlogit, n_snps)
    fit = rfe.fit(X_std, y)
    #sub_preds = preds[fit.support_]
    pred_str = ' + '.join(preds[fit.support_])
    X_sub = X_std[:, fit.support_]

    for j in range(n_reps):
        k_fold = KFold(n_splits = val_k, random_state = j, shuffle = True)
        sub_score_arr[j,:] = cross_val_score(mlogit, X_sub, y, cv = k_fold, scoring = 'accuracy')

    filled_df = proto_df
    filled_df['accuracy'] = sub_score_arr.flatten()
    filled_df['nSNPs'] = n_snps
    filled_df['SNP_IDs'] = pred_str



