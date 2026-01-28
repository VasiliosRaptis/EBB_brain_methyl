#!/bin/bash

## initialise the modules framework
source /etc/profile.d/modules.sh

## load modules
module load igmm/apps/plink/1.90b7.2
module load igmm/apps/plink/2.00a6LM

## variables
PLINK_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/plink_files
LOGS_DIR=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/logs
MAPPING_FILE=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/meffil_data/sampleID.mapping.genconc.txt

## update IDs
cat $MAPPING_FILE | awk 'NR>1 && $3 >= 0.75 {print $1,$1,$2,$2}' > $PLINK_DIR/update_ids.txt
cat $MAPPING_FILE | awk 'NR>1 && $3 >= 0.75 {print $2,$2}' >  $PLINK_DIR/good_samples.txt

plink2 \
  --pfile $PLINK_DIR/pgen/imputed_allchr \
  --update-ids $PLINK_DIR/update_ids.txt \
  --keep $PLINK_DIR/good_samples.txt \
  --make-pgen \
  --out $PLINK_DIR/pgen/imputed_allchr_newIDs \
  --threads $(nproc)

## clean-up
mv $PLINK_DIR/pgen/*.log $LOGS_DIR