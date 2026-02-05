#!/usr/bin/env python3


import os
import glob
import pandas as pd
from pandas.core.common import flatten
import sys

batch = sys.argv[1]
print(batch)

# list of all CpG names
sspaths = os.listdir('scores/scores/')
cpgs = list({x.replace('.sscore', '').replace('.tmp', '') for x in sspaths})
n_cpgs = len(cpgs)
print(n_cpgs)
cpgs[0:4]

# initiate predict.txt matrix: samples x scores
input0 = 'scores/scores/' + cpgs[0] + '.sscore'
pred_mat = pd.read_csv(input0, sep="\t").iloc[:,[0,1]]
pred_mat.columns = ['FID', 'IID']
pred_mat.set_index('IID', drop=True, inplace=True)

for i in range(n_cpgs):
    # get cpg i
    cpg_i = cpgs[i]
    # get score for cpg i
    input_i = 'scores/scores/' + cpg_i + '.sscore'
    try:
        score_i = pd.read_csv(input_i, sep="\t").iloc[:,[0,4]]
    except Exception as e:
        print(f"Skipping {cpg_i} due to error: {e}")
        continue
    score_i.columns = ['IID', cpg_i]
    score_i.set_index('IID', drop=True, inplace=True)
    # add to matrix 
    pred_mat = pred_mat.join(score_i, on='IID', how='inner')
    
pred_mat.reset_index(inplace=True)
print(pred_mat.shape)

# save
out = 'batch' + batch + '_predict.txt'
pred_mat.to_csv(out, header=True, index=False, sep='\t')