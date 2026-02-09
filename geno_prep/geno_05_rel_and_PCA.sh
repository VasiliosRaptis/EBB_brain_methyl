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



### get PCA outliers & # of PC to include
Rscript geno_05.1_rel_and_PCA_helper.R

## exclude outlier
plink2  \
  --pfile $PLINK_DIR/pgen/imputed_allchr_newIDs \
  --remove $PLINK_DIR/rel_mat/outliers.txt \
  --make-pgen  \
  --out $PLINK_DIR/pgen/imputed_allchr_newIDs_noOutlier  \
  --threads $(nproc)
