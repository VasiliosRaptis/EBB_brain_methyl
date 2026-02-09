## load libraries
library(stringr)
library(data.table) 
library(vroom)
library(ggplot2)
#library(tidyr)
#library(limma)
#library(meffil)
library(readxl)
library(dplyr)

# set wd
setwd('/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/')

## load files
# normalised methulation betas (phenotypes)
load('meffil_data/norm.beta.pc10clean.Robj')

# make new matrix
norm.beta.int = matrix(nrow = nrow(norm.beta), ncol = ncol(norm.beta))
rownames(norm.beta.int) <- rownames(norm.beta)
colnames(norm.beta.int) <- colnames(norm.beta)

for (i in 1:nrow(norm.beta)) {
    cpg_i = rownames(norm.beta)[i]
    # get phenotype for cpg_i
    norm.beta_i = norm.beta[cpg_i,] 
    # int transform cpg_i
    int_i = qnorm((rank(norm.beta_i,na.last="keep")-0.5)/sum(!is.na(norm.beta_i)))
    norm.beta.int[cpg_i,] <- int_i
}

## save
save(norm.beta.int, file='meffil_data/norm.beta.pc10clean.int.Robj')