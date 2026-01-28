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
module load igmm/apps/bcftools/1.6
module load igmm/apps/plink/2.00a6LM

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

### notes: 
### - requires to have downloaded and unzipped per chromosome imputed vcf files from michigan imputation server 
### - bcftools: filter based on INFO score  
### - plink2: will convert vcf -> .pgen (plink2 better at handling vcf files (see: https://groups.google.com/g/plink2-users/c/GxBcwK3Qw6w) & do QC 
### - plink2: merge all chrom



INFO=0.7
mkdir -p $PLINK_DIR/pgen
mkdir -p $PLINK_DIR/extract

# change the update sex file beacuse the FIDs/IIDs changed (will change back later)
cat $PLINK_DIR/update_sex.txt | awk '{print "0_"$1,"0_"$2,$3}' > $PLINK_DIR/update_sex_forpgen.txt 

for CHROM in {1..22}  
do
  ## extract list of typed variants and with INFO threshold
  bcftools filter -i "INFO/TYPED==1 || INFO/R2 >= $INFO" $MIS_DIR/results/chr${CHROM}.info.gz |  bcftools query -f '%ID\n' > $PLINK_DIR/extract/chr${CHROM}_extract.txt
  ## convert to plink .pgen & QC 
  plink2 \
    --vcf $MIS_DIR/results/chr${CHROM}.dose.vcf.gz \
    --double-id \
    --extract $PLINK_DIR/extract/chr${CHROM}_extract.txt \
    --maf 0.01 \
    --mac 3 \
    --geno 0.05 --hwe 1e-6 --mind 0.05 \
    --make-pgen --out $PLINK_DIR/pgen/chr${CHROM} \
    --update-sex $PLINK_DIR/update_sex_forpgen.txt \
    --threads $(nproc) 

  echo "Finished $CHROM"
  
done 

for CHROM in X 
do
  ## extract list of typed variants and with INFO threshold
  bcftools filter -i "INFO/TYPED==1 || INFO/R2 >= $INFO" $MIS_DIR/results/chr${CHROM}.info.gz |  bcftools query -f '%ID\n' > $PLINK_DIR/extract/chr${CHROM}_extract.txt
  ## convert to plink .pgen & QC 
  plink2 \
    --vcf $MIS_DIR/results/chr${CHROM}.dose.vcf.gz \
    --double-id \
    --extract $PLINK_DIR/extract/chr${CHROM}_extract.txt \
    --maf 0.01 \
    --mac 3 \
    --geno 0.05 --hwe 1e-6 --mind 0.05 \
    --make-pgen --out $PLINK_DIR/pgen/chr${CHROM} \
    --update-sex $PLINK_DIR/update_sex_forpgen.txt \
    --split-par hg38\
    --threads $(nproc) 

  echo "Finished $CHROM"
  
done 

## merge .pgen files - set IDs to CHROM:POS:REF:ALT

ls -1 $PLINK_DIR/pgen/chr*.pgen | sed 's/\.pgen$//' > $PLINK_DIR/pgen/pgen-merge-list.txt
cat $PLINK_DIR/pgen/chr22.psam | awk 'NR>1 {$3 = gensub(/^0_/, "", 1, $1); print $1,$2,$3,$3}' > $PLINK_DIR/update_ids_forpgen.txt

### NOTE: run with plink2 version in ../../software/plink2/plink2, otherwise threse is a bug in merging 
plink2 \
#../../software/plink2/plink2 \
  --pmerge-list $PLINK_DIR/pgen/pgen-merge-list.txt \
  --set-all-var-ids '@:#:$r:$a' \
  --new-id-max-allele-len 342 \
  --update-ids $PLINK_DIR/update_ids_forpgen.txt \
  --make-pgen \
  --out $PLINK_DIR/pgen/imputed_allchr

## clean-up
rm $PLINK_DIR/pgen/imputed_allchr-merge*
rm $PLINK_DIR/pgen/chr*
mv $PLINK_DIR/pgen/*.log $LOGS_DIR








