#!/usr/bin/env python3

import os
import glob
import pandas as pd
import sys
import statsmodels.api as sm
import numpy as np
from scipy import stats
from scipy.stats import zscore
from mwas_func import get_pheno_df

## passed arguments
batch   = sys.argv[1] # 1,2,3 .. 50 -> from for loop in .sh 
phe_map = sys.argv[2] # phecode mapping file
icdfile = sys.argv[3] # cases file
covfile = sys.argv[4] # covariates/controls file
outsfx  = sys.argv[5] # output df filename suffix

## load predictions files - for all phenotypes
print('loading predictions')
predfile = 'batch' + str(batch) + '_predict.txt'
pred = pd.read_csv(predfile, sep="\t").rename(columns={'IID':'eid'}).drop(columns=['FID']).set_index('eid')

## load phecodes mapping file
phecodes = pd.read_csv(phe_map)
phecodes = phecodes[['phecode', 'description']].drop_duplicates()
n_pheno  = phecodes.shape[0] 

## set covariates
PCS = [f"PC{i}" for i in range(1, 21)]
covariates = ['sex', 'age', 'batch', 'e4_status'] + PCS

for j in range(n_pheno):
    
    ## start phenotype j analysis
    pheno_name_j  = phecodes.iloc[j]['phecode']
    description_j = phecodes.iloc[j]['description']
    
    ## load phenotype j file
    print(f'starting analysis for {description_j} ({pheno_name_j}) | {j} / {n_pheno}:')
    print(f'\tloading phenotype/covariates')

    phe_j = get_pheno_df(
        pheno=pheno_name_j,
        icdfile=icdfile,
        covfile=covfile,
        covariates = covariates
    )
    
    ## match the order of pred matrix eids
    phe_j = phe_j.reindex(pred.index)

    ## Get phenotype residuals
    print('\tResidualising phenotype')

    Y_orig  = phe_j[pheno_name_j]
    X_covar = sm.add_constant(phe_j[covariates])
    Y_resid = sm.Logit(Y_orig, X_covar).fit(disp=0).resid_response

    ## get n cases /controls
    n_cas = len(Y_orig[Y_orig==1])
    n_cnt = len(Y_orig[Y_orig==0])

    ## Association analysis
    cpgs   = pred.columns[2:].tolist()
    n_cpgs = len(cpgs)
    print(f"\t{description_j} ~ CpGs associations; {n_cpgs} CpGs")

    res = []

    y = Y_resid.values

    ## pheno i resid ~ CpG
    for i in range(n_cpgs):

        # get scaled CpG i data
        cpg_i = cpgs[i]  
        x = zscore(pred[cpg_i], nan_policy='omit')

        # drop missing
        mask = ~np.isnan(x) & ~np.isnan(y)
        x_ = x[mask]
        y_ = y[mask]
        n = len(x_)

        if n < 3:
            continue

        # Fit model
        slope, intercept, r, pval, stderr = stats.linregress(x_, y_)

        # 95% CI using t distribution
        tcrit = stats.t.ppf(0.975, df=n - 2)
        lci = slope - tcrit * stderr
        uci = slope + tcrit * stderr

        res.append({
            "cpg": cpg_i,
            "beta": slope,
            "se": stderr,
            "lci": lci,
            "uci": uci,
            "zscore": slope / stderr,
            "pvalue": pval,
            "r": r,
            "n_cas": n_cas,
            "n_cnt": n_cnt
        })

    res_df = pd.DataFrame(res)

    ## save results df
    print('\tSaving results file')
    outfile_j = pheno_name_j + '_' + outsfx
    res_df.to_csv(outfile_j, index=False, sep = ",")
