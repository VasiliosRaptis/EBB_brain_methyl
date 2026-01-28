#!/bin/bash

## initialise the modules framework
source /etc/profile.d/modules.sh

## load modules
module load igmm/apps/plink/1.90b7.2
module load igmm/apps/plink/2.00a6LM

## variables
PLINK_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/plink_files
LOGS_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/logs

## note: see https://www.cog-genomics.org/plink/2.0/tutorials/qc2b

## get non imputed set 
cat $PLINK_DIR/pgen/imputed_allchr_newIDs.pvar | grep 'TYPED' | awk 'NR >1 {print $3}' > $PLINK_DIR/non_imputed.txt

plink2  \
  --pfile $PLINK_DIR/pgen/imputed_allchr_newIDs \
  --extract $PLINK_DIR/non_imputed.txt  \
  --make-pgen  \
  --out $PLINK_DIR/pgen/non_imputed_allchr_newIDs  \
  --threads $(nproc)

## ld-pruning
mkdir $PLINK_DIR/rel_mat
r2=0.2

plink2 \
  --pfile $PLINK_DIR/pgen/non_imputed_allchr_newIDs  \
  --autosome \
  --snps-only \
  --indep-pairwise 500kb $r2 \
  --out $PLINK_DIR/rel_mat/non_imputed_postqc

plink2 \
  --pfile $PLINK_DIR/pgen/non_imputed_allchr_newIDs \
  --make-pgen \
  --extract $PLINK_DIR/rel_mat/non_imputed_postqc.prune.in  \
  --out $PLINK_DIR/rel_mat/non_imputed_postqc_pruned${r2}

## make relationship matrix

plink2 \
  --pfile $PLINK_DIR/rel_mat/non_imputed_postqc_pruned${r2}  \
  --make-king square0 \
  --maf 0.1 \
  --out $PLINK_DIR/rel_mat/non_imputed_postqc


## PCA
plink2 \
  --pfile $PLINK_DIR/rel_mat/non_imputed_postqc_pruned${r2} \
  --freq counts \
  --keep-founders \
  --pca biallelic-var-wts \
  --out $PLINK_DIR/rel_mat/non_imputed_postqc_pca



### in R:
setwd("/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata")
library(dplyr)
library(data.table)
library(ggplot2)

## load king matrix and flatten 
king <- fread('plink_files/rel_mat/non_imputed_postqc.king') %>% as.matrix
id <- fread('plink_files/rel_mat/non_imputed_postqc.king.id')
colnames(king) <- id$IID
rownames(king) <- id$IID

flattenKingMatrix <- function(cormat) {
  lt <- lower.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[lt]],
    column = rownames(cormat)[col(cormat)[lt]],
    coef  = (cormat)[lt]
    )
}

king_df <- flattenKingMatrix(king)

## summary (no relatives)
summary(king_df$coef)

## plot king coefficients 
png('plink_files/rel_mat/king_plot.png')
hist(king_df$coef, xlim=c(-.05,0.5))
dev.off()

## how many PCs to include?
eigenval <- fread('plink_files/rel_mat/non_imputed_postqc_pca.eigenval')

## make scree plot -> eigenvalues flatten at PC=4
png('plink_files/rel_mat/eigenval_plot.png')
plot(eigenval$V1, main='GRM Scree plot')
dev.off()



## find population outliers (+-5df from POP PC1&2 mean)
pca <- fread('plink_files/rel_mat/non_imputed_postqc_pca.eigenvec')

pc1_min <- mean(pca$PC1) - sd(pca$PC1)*5
pc1_max <- mean(pca$PC1) + sd(pca$PC1)*5
pc2_min <- mean(pca$PC2) - sd(pca$PC2)*5
pc2_max <- mean(pca$PC2) + sd(pca$PC2)*5

pca2 <-
pca %>% 
  mutate(outlier = ifelse( PC1 < pc1_min | PC1 > pc1_max | PC2 < pc2_min | PC2 > pc2_max, 1, 0)) %>% 
  mutate(outlier = as.factor(outlier))

png('plink_files/rel_mat/pca_plot.png')
ggplot(pca2, aes(x=PC1, y=PC2, color=outlier)) + geom_point() + theme_bw()
dev.off()

## "SD037/14B" is population outlier to remove
outlier.ids <- pca2 %>% filter(outlier == 1) %>% pull(IID) 


## exclude outlier
plink2  \
  --pfile $PLINK_DIR/pgen/imputed_allchr_newIDs \
  --remove $PLINK_DIR/rel_mat/outliers.txt \
  --make-pgen  \
  --out $PLINK_DIR/pgen/imputed_allchr_newIDs_noOutlier  \
  --threads $(nproc)
