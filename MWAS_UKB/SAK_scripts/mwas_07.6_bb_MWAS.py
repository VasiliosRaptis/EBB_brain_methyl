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
print('loading brain-based predictions ...')
brain_predfile = 'batch' + str(batch) + '_predict.txt'
brain_pred = pd.read_csv(brain_predfile, sep="\t").rename(columns={'IID':'eid'}).drop(columns=['FID']).set_index('eid')

print('loading blood-based predictions ...')
blood_predfile = 'batch' + str(batch) + '_predict_blood.txt'
blood_pred = pd.read_csv(blood_predfile, sep="\t").rename(columns={'IID':'eid'}).drop(columns=['FID']).set_index('eid')

## get common cpgs 
common_cpgs = list(blood_pred.columns.intersection(brain_pred.columns))
n_common_cpgs = len(common_cpgs)
print(f'{n_common_cpgs} CpG sites both in brain- & blood-based DNAm predictors')

## load phecodes mapping file
phecodes = pd.read_csv(phe_map)
phecodes = phecodes[['phecode', 'description']].drop_duplicates()
n_pheno  = phecodes.shape[0] 

## set covariates
PCS = [f"PC{i}" for i in range(1, 21)]
covariates = ['sex', 'age', 'batch'] + PCS

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
    phe_j = phe_j.reindex(brain_pred.index)
    blood_pred = blood_pred.reindex(brain_pred.index)

    ## Get phenotype residuals
    print('\tResidualising phenotype')

    Y_orig  = phe_j[pheno_name_j]
    X_covar = sm.add_constant(phe_j[covariates])
    Y_resid = sm.Logit(Y_orig, X_covar).fit(disp=0).resid_response

    ## get n cases /controls
    n_cas = len(Y_orig[Y_orig==1])
    n_cnt = len(Y_orig[Y_orig==0])

    ## Association analysis
    cpgs   = list(blood_pred.columns.intersection(brain_pred.columns))
    n_cpgs = len(cpgs)
    print(f"\t{description_j} ~ CpGs associations; {n_cpgs} CpGs")

    res = []

    y = Y_resid.values

    ## pheno i resid ~ CpG
    for i in range(n_cpgs):

        ###  get scaled CpG i data
        cpg_i = cpgs[i]

        # blood-based predictor
        x1 = zscore(blood_pred[cpg_i], nan_policy='omit')
        # x1 = blood_pred[cpg_i]

        # brain-based predictor 
        x2 = zscore(brain_pred[cpg_i], nan_policy='omit')
        # x2 = brain_pred[cpg_i]

        # mask missing
        mask = ~np.isnan(x1) & ~np.isnan(x2) & ~np.isnan(y)
        x1_, x2_, y_ = x1[mask], x2[mask], y[mask]
        n = len(y_)

        if n < 3:
            continue

        # design matrix
        X = np.column_stack((x1_, x2_))
        X = sm.add_constant(X)

        model = sm.OLS(y_, X).fit()

        # coefficients
        beta1, beta2 = model.params[1], model.params[2]
        se1, se2     = model.bse[1], model.bse[2]
        pval1, pval2 = model.pvalues[1], model.pvalues[2]

        ci = model.conf_int()

        # partial correlations 
        sd_y  = np.std(y_, ddof=1)
        sd_x1 = np.std(x1_, ddof=1)
        sd_x2 = np.std(x2_, ddof=1)

        r_x1 = beta1 * sd_x1 / sd_y
        r_x2 = beta2 * sd_x2 / sd_y

        # brain-blood predictor correlation 
        r_x12 = stats.pearsonr(x1_, x2_)[0]

        res.append({
            "cpg": cpg_i,

            # x1
            "beta_x1": beta1,
            "se_x1": se1,
            "lci_x1": ci[1, 0],
            "uci_x1": ci[1, 1],
            "zscore_x1": beta1 / se1,
            "pvalue_x1": pval1,
            "r_x1": r_x1,

            # x2
            "beta_x2": beta2,
            "se_x2": se2,
            "lci_x2": ci[2, 0],
            "uci_x2": ci[2, 1],
            "zscore_x2": beta2 / se2,
            "pvalue_x2": pval2,
            "r_x2": r_x2,

            "n_cas": n_cas,
            "n_cnt": n_cnt,

            "r_x12": r_x12

        })
    
    res_df = pd.DataFrame(res)

    ## save results df
    print('\tSaving results file')
    outfile_j = pheno_name_j + '_' + outsfx
    res_df.to_csv(outfile_j, index=False, sep = ",")
