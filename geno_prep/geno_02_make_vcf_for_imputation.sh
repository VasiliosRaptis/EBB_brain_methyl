#!/bin/bash

## initialise the modules framework
source /etc/profile.d/modules.sh

## load modules
module load igmm/apps/plink/1.90b7.2
module load igmm/apps/bcftools/1.6

## variables
PLINK_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/plink_files
VCF_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/vcf_files
LOGS_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/logs
MIS_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/mis_outputs

## setup
mkdir -p $PLINK_DIR
mkdir -p $VCF_DIR
mkdir -p $LOGS_DIR
mkdir -p $MIS_DIR

### note: for michigan imputation server input https://imputationserver.sph.umich.edu/index.html#!

## make vcf
plink \
  --bfile $PLINK_DIR/methylation_postqc \
  --recode vcf-iid \
  --out $VCF_DIR/methylation_postqc \
  --threads $(nproc)

## make vcf.gz + index file
bgzip -c $VCF_DIR/methylation_postqc.vcf > $VCF_DIR/methylation_postqc.vcf.gz
tabix -p vcf $VCF_DIR/methylation_postqc.vcf.gz

## split by chromosome (see https://www.biostars.org/p/9506417/#9555305)
bcftools index -s $VCF_DIR/methylation_postqc.vcf.gz | cut -f 1 | while read C; do bcftools view -O z -o $VCF_DIR/split.${C}.vcf.gz $VCF_DIR/methylation_postqc.vcf.gz "${C}" ; done

## switch incorrect alleles (based on MIS error output) + change X chr code & and create vcf files again
cat $MIS_DIR/reports/snps-excluded*.txt | awk -F'\t' '$6 ~ /Allele switch/ {print $1,$4,$5}' > $PLINK_DIR/a2_switch.txt

plink \
  --vcf $VCF_DIR/methylation_postqc.vcf \
  --const-fid \
  --a2-allele $PLINK_DIR/a2_switch.txt 3 1 \
  --output-chr M \
  --recode vcf \
  --out $VCF_DIR/methylation_postqc_switch

bgzip -c $VCF_DIR/methylation_postqc_switch.vcf > $VCF_DIR/methylation_postqc_switch.vcf.gz
tabix -p vcf $VCF_DIR/methylation_postqc_switch.vcf.gz

rm $VCF_DIR/split.*
bcftools index -s $VCF_DIR/methylation_postqc_switch.vcf.gz | cut -f 1 | while read C; do bcftools view -O z -o $VCF_DIR/split.${C}.vcf.gz $VCF_DIR/methylation_postqc_switch.vcf.gz "${C}" ; done

## clean-up
mv $PLINK_DIR/*.log $LOGS_DIR
mv $VCF_DIR/*.log $LOGS_DIR
