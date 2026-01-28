#!/bin/bash

#$ -cwd
#$ -N score
#$ -o ./qsub_logs/score.$JOB_ID.$TASK_ID.o
#$ -e ./qsub_logs/score.$JOB_ID.$TASK_ID.e
#$ -l h_rt=47:00:00
#$ -l h_rss=5G
##$ -pe sharedmem 4
#$ -t 1-50
#$ -tc 50

### initialise the modules framework
source /etc/profile.d/modules.sh

### setup chunk number
if (( $SGE_TASK_ID >= 1 && $SGE_TASK_ID <= 9 )); then
    chunk=0$SGE_TASK_ID
else
    chunk=$SGE_TASK_ID
fi

### Setup
module load igmm/apps/plink/1.90b7.2
module load igmm/apps/plink/2.00a6LM

software=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/software
EBB=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation
WD=$EBB/MWAS
INP=$WD/INP
TMP=$WD/TMP
OUT=$WD/SCORES

# inputs
PHENO=$EBB/tensorqtl/phenotypes.int.bed
EBB_bfile=$INP/ldref
WGTLIST_FULL=$WD/EBB.BRAIN.METHYL.HERIT.list # list of paths to weights (hsq < 0.05)

cd $WD
mkdir -p $OUT/input
mkdir -p $OUT/scores/EBB/summary
mkdir -p $OUT/scores/EBB/logs

### split weights list into chunks - RUN ONCE
# split -n l/50 --numeric-suffixes=1 $WGTLIST_FULL $TMP/temp_WGTLIST_

### BATCH STARTS HERE

### make imputed methyltation scores for each CpG 
### use the *score.in files to impute methylation in UKB 
### see: https://github.com/gusevlab/fusion_twas/tree/master/utils

WGTLIST=$TMP/temp_WGTLIST_$chunk

a=($(wc -l $WGTLIST))
n_cpgs=${a[0]} # total no. of cpgs 

# initiate pheno ~ score correlation table
CORRTABLE=$OUT/scores/EBB/summary/scores.cor.table.chunk$chunk
>$CORRTABLE

echo "Analysis started at: $(date)"

# set timer
total_time=0
n_loops=0

for i in $(seq 1 "$n_cpgs"); do
  start=$(date +%s.%N)
  
  # pull path to weights
  wgtfile_i=$(awk -v i=$i 'NR==i' $WGTLIST)
  # pull cpg name
  cpg_name_i=$(basename $wgtfile_i .wgt.RDat)
  # pull phenotype for cpg i and make pheno file for plink
  awk -F'\t' -v cpg=$cpg_name_i '$4==cpg{print; exit}' $PHENO | tr '\t' '\n' | tail -n+5 | paste -d" " $INP/samples.id - > $TMP/$cpg_name_i.pheno.for.score
  # make per-cpg input for plink2 --score
  Rscript $software/fusion_twas-master/utils/make_score.R $wgtfile_i 1> $OUT/input/$cpg_name_i.score.in 
  # make per-cpg score (in EBB) + add phenotype
  plink2 --bfile $EBB_bfile --pheno $TMP/$cpg_name_i.pheno.for.score --score $OUT/input/$cpg_name_i.score.in 1 2 4 header-read --out $OUT/scores/EBB/$cpg_name_i --threads $(nproc) &>/dev/null
  # get pheno ~ score correlation (in EBB)
  Rscript $EBB/scripts/x_MWAS_04.3b_make_scores_helper.R $OUT/scores/EBB/$cpg_name_i.sscore $cpg_name_i >> $CORRTABLE
  # cleanup
  rm $TMP/$cpg_name_i.pheno.for.score
  mv $OUT/scores/EBB/$cpg_name_i.log $OUT/scores/EBB/logs/
  
  # print progress
  pct=$(echo "100 * $i / $n_cpgs" | bc -l)
  
  end=$(date +%s.%N)
  time_i=$(echo "$end - $start" | bc -l)
  
  n_loops=$((n_loops + 1))
  total_time=$(echo "$total_time + $time_i" | bc -l)
  avg_time=$(echo "$total_time / $n_loops" | bc -l)
  
  printf "finished CpG %s, %d out of %d (%.2f%%), in %.2f s | avg time per CpG: %.2f s\n" "$cpg_name_i" "$i" "$n_cpgs" "$pct" "$time_i" "$avg_time"
  
done

echo "Analysis finished at: $(date)"

### merge all correlation tables
chunks_finished=$(ls -1 $OUT/scores/EBB/summary/scores.cor.table.chunk* | wc -l)
if (( $chunks_finished == 49 )); then
  echo -e "cpg\tmodel\tcorr\tcorr.pv\tcorr.lci\tcorr.uci\tr2" > $OUT/scores/EBB/summary/scores.cor.table.txt
  cat $OUT/scores/EBB/summary/scores.cor.table.txt $OUT/scores/EBB/summary/scores.cor.table.chunk* >> $OUT/scores/EBB/summary/scores.cor.table.txt
fi

