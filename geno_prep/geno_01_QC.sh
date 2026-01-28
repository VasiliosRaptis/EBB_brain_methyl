#!/bin/bash

#$ -cwd
#$ -N EBB_geno
#$ -o ./qsub_logs/QC.$TASK_ID.o
#$ -e ./qsub_logs/QC.$JOB_ID.$TASK_ID.e
#$ -l h_rt=17:00:00
#$ -l h_vmem=4G
#$ -pe sharedmem 4
##$ -t 1-22
##$ -tc 22

## initialise the modules framework
source /etc/profile.d/modules.sh

## load modules
module load igmm/apps/plink/1.90b7.2
#module load roslin/regenie/3.2.2

## variables
RAW_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/raw_files
OUT_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/plink_files
LOGS_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/logs


## setup
mkdir -p $OUT_DIR
mkdir -p $LOGS_DIR

## convert ped/map to bed/bim/fam
#plink \
#  --file $RAW_DIR/methylation \
#  --keep-allele-order \
#  --no-sex --no-parents --no-pheno --no-fid \
#  --make-bed \
#  --out $OUT_DIR/methylation \
#  --threads $(nproc) 

## allele frequency report
#plink --bfile $OUT_DIR/methylation --freq --out $OUT_DIR/methylation

## update sex
plink  \
  --bfile $OUT_DIR/methylation \
  --check-sex \
  --make-bed \
  --out $OUT_DIR/methylation_sex \
  --update-sex $OUT_DIR/update_sex.txt \
  --threads $(nproc) 

## QC
plink \
  --bfile $OUT_DIR/methylation_sex \
  --maf 0.01 \
  --mac 3 \
  --geno 0.05 --hwe 1e-6 --mind 0.05 \
  --make-bed \
  --out $OUT_DIR/methylation_postqc \
  --threads $(nproc)

## clean-up
mv $OUT_DIR/*.log $LOGS_DIR