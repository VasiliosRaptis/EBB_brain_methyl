#!/usr/bin/env python3

import pandas as pd
import numpy as np

def get_pheno_df(pheno, icdfile, covfile, covariates=None):
    """
    Generate a case-control phenotype DataFrame for a specific phecode.

    Parameters
    ----------
    pheno : str
        The phecode of interest.
    icdfile : str
        Path to the ICD file (must have 'eid', 'phecode', 'date').
    covfile : str
        Path to the covariate file (must have 'eid', 'dob', 'sex', 'batch', 'e4_status', PC1-PC20).
    covariates : list of str, optional
        Covariates to include. Default is ['sex', 'age', 'batch', 'PC1' to 'PC20'].

    Returns
    -------
    pheno_i : pd.DataFrame
        DataFrame with cases (phecode==1) and controls (phecode==0) with selected covariates.
    """
    
    # Default covariates
    if covariates is None:
        PCS = [f"PC{i}" for i in range(1, 21)]
        covariates = ['sex', 'age', 'batch'] + PCS
    
    # Load data
    icd_phe = pd.read_csv(icdfile).set_index('eid')
    cov     = pd.read_csv(covfile).set_index('eid')
    
    # Get cases
    icd_i = icd_phe.loc[icd_phe['phecode'] == pheno]
    icd_i_eids = np.unique(icd_i.index.values)
    
    cas_i = cov.loc[cov.index.intersection(icd_i_eids)].copy()
    
    # keep earliest diagnosis date 
    icd_i = icd_i.reset_index().sort_values(by=['eid', 'date']).set_index('eid')
    icd_i = icd_i[ ~ icd_i.index.duplicated(keep='first') ]

    # Add earliest diagnosis date and compute age
    cas_i = cas_i.join(icd_i['date'], how='inner')
    cas_i['date'] = pd.to_datetime(cas_i['date'])
    cas_i['dob'] = pd.to_datetime(cas_i['dob'])
    cas_i['age'] = ((cas_i['date'] - cas_i['dob']).dt.days / 365.25).round(1)
    
    # Keep selected covariates
    cas_i = cas_i.loc[:, covariates]
    
    # Add phecode label
    cas_i[pheno] = 1
    
    # Get controls
    clt_i = cov.loc[~cov.index.isin(icd_i_eids)].copy()
    
    # Rename censoring age if present
    if 'censoring_age' in clt_i.columns and 'age' in covariates:
        clt_i = clt_i.rename(columns={'censoring_age': 'age'})
    
    # Keep selected covariates
    clt_i = clt_i.loc[:, covariates]
    
    # Add phecode label
    clt_i[pheno] = 0
    
    # Combine cases and controls
    pheno_i = pd.concat([cas_i, clt_i], ignore_index=False)
    
    return pheno_i