library(qvalue)
library(dplyr)
library(data.table)

## adapted from: https://github.com/broadinstitute/tensorqtl/blob/0c4db65a0cdc47f3b824ae530b89d270ef5e0096/tensorqtl/post.py#L18

# load cis results
cis <- fread('tensorqtl/results/EBB.cis_qtl_perm_all.csv')

print('Computing q-values:')
print(paste('Number of phenotypes tested:', nrow(cis)))

# check p-values
r =  cor(cis$pval_beta, cis$pval_perm, method='p')
print(paste('Correlation between Beta-approximated and empirical p-values:', round(r,5)))


## calculate qvalues
fdr = 0.05

qobj    <- qvalue(cis$pval_beta, fdr=fdr)
qvalues <- qobj$qvalues
pi0     <- qobj$pi0
lfdr    <- qobj$lfdr

cis$qval <- qvalues

print(paste('Proportion of significant phenotypes (1-pi0):', 1-pi0))
print(paste('QTL phenotypes @ FDR', fdr, ':', nrow(cis[cis$qval<fdr,])))

# determine global min(p) significance threshold and calculate nominal p-value threshold for each gene


lb <- sort(cis[cis$qval <= fdr,]$pval_beta)
ub <- sort(cis[cis$qval >  fdr,]$pval_beta)

if (length(lb) > 0) {
  lb = lb[length(lb)]
  if (length(ub) > 0) {
    ub = ub[1]
    pthreshold = (lb+ub)/2
  }
  else pthreshold = lb
}

print(paste('min p-value threshold @ FDR', fdr, pthreshold))

cis$pval_nominal_threshold <- qbeta(pthreshold, cis$beta_shape1, cis$beta_shape2)

# save
fwrite(cis, 'tensorqtl/results/EBB.cis_qtl_perm_all_qval.csv')
