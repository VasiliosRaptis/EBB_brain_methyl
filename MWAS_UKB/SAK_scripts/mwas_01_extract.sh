#!/bin/bash

# $1 chr -> passed from dx run swiss-army-knife
CHR=$1
KEEP=wb.eids
EXTRACT=imp_c${CHR}.extract
OUT=imp_wb_qc_c${CHR}

plink2 \
  --bgen ukb22828_c${CHR}_b0_v3.bgen ref-first \
  --sample ukb22828_c${CHR}_b0_v3.sample \
  --maf 0.01 \
  --mac 3 \
  --geno 0.05 --hwe 1e-10 --mind 0.05 \
  --keep $KEEP \
  --extract $EXTRACT \
  --make-bed \
  --out ${OUT} \
  --threads $(nproc)
