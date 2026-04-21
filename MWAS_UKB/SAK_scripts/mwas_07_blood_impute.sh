#!/bin/bash

### inputs needed: 
## passed by -iin: EPICDB.BLOOD.METHYL.HERIT.tar.bz2, imp_wb_qc_newid_all* plink files 
## passed by -iin: EBB.BRAIN.METHYL.HERIT.list
## passed by -icmd: $1 batch number

### setup batch number
# $1 batch -> passed from dx run swiss-army-knife
batch=$1

### setup

# extract weights
tar -xjf EPICDB.BLOOD.METHYL.HERIT.tar.bz2 2>/dev/null

# inputs
UKB_bfile=imp_wb_qc_newid_all
WGTLIST_FULL_EPICDB=EPICDB.BLOOD.METHYL.HERIT/scores/scores.list # list of paths to weights (hsq < 0.05)
WGTLIST_FULL_EBB=EBB.BRAIN.METHYL.HERIT.list
scoresdir=EPICDB.BLOOD.METHYL.HERIT/scores/input

# output
mkdir -p scores/scores
mkdir -p scores/logs
>errors_$batch.log # capture errors 
>times_$batch.log # print times

### BATCH SETUP STARTS HERE

## get same CpGs per batch as in the brain-based analysis 

# ebb 
mkdir -p TMP_EBB
split -n l/50 --numeric-suffixes=1 $WGTLIST_FULL_EBB TMP_EBB/temp_WGTLIST_
for f in TMP_EBB/temp_WGTLIST_0*; do
  mv $f $(echo $f | sed 's/_0/_/')
done
WGTLIST_EBB=TMP_EBB/temp_WGTLIST_$batch # brain-based CpG list

codes=$(sed 's|.*/\(.*\)\.wgt\.RDat|\1|' $WGTLIST_EBB)
pattern=$(echo "$codes" | tr '\n' '|')
pattern="${pattern%|}"  # Remove the trailing pipe

# epic-db
mkdir -p TMP_EPICDB
grep -E "$pattern" $WGTLIST_FULL_EPICDB > TMP_EPICDB/temp_WGTLIST_$batch
WGTLIST=TMP_EPICDB/temp_WGTLIST_$batch # to be passed by SAK

### BATCH SETUP ENDS HERE


### make imputed methyltation scores for each CpG 
### use the *score.in files to impute methylation in UKB 
### see: https://github.com/gusevlab/fusion_twas/tree/master/utils

#a=($(wc -l $WGTLIST))
#n_cpgs=${a[0]} # total no. of cpgs 
n_cpgs=$(wc -l $WGTLIST | awk '{print $1}') # total no. of cpgs 

echo "Analysis started at: $(date)"

# set timer
total_time=0
n_loops=0

for i in $(seq 1 "$n_cpgs"); do
  start=$(date +%s.%N)
  
  # pull path to weights
  wgtfile_i=$(awk -v i=$i 'NR==i' $WGTLIST)
  # pull cpg name
  cpg_name_i=$(basename $wgtfile_i .score.in)
  # make per-cpg score in UKB; note: use --extract otherwise allele freq calcuation takes for ever
  plink2 \
    --bfile $UKB_bfile \
    --extract $scoresdir/$cpg_name_i.score.in \
    --score $scoresdir/$cpg_name_i.score.in 1 2 4 header-read \
    --out scores/scores/$cpg_name_i \
    --threads $(nproc)  >/dev/null 2>> errors_$batch.log
  # cleanup
  mv scores/scores/$cpg_name_i.log scores/logs
  
  # print progress
  pct=$(awk "BEGIN { printf \"%.2f\", 100 * $i / $n_cpgs }")

  end=$(date +%s.%N)
  time_i=$(awk "BEGIN { printf \"%.2f\", $end - $start }")

  n_loops=$(( n_loops + 1 ))
  total_time=$(awk "BEGIN { print $total_time + $time_i }")
  avg_time=$(awk "BEGIN { printf \"%.2f\", $total_time / $n_loops }")

  printf "finished CpG %s, %d out of %d (%s%%), in %.2f s | avg time per CpG: %.2f s\n" \
         "$cpg_name_i" "$i" "$n_cpgs" "$pct" "$time_i" "$avg_time"  >> times_$batch.log
done

echo "Analysis finished at: $(date)"

# make matrix of imputed methyltation scores: samples x scores 
echo "samples x scores matrix started: $(date)"
pip install pandas
python3 mwas_02.1_pred_mat.py $batch
echo "samples x scores matrix finished: $(date)"

# make .tar.bz2 file 
echo "compressing started at: $(date)"
tar -cjSf scores_$batch.tar.bz2 scores
echo "compressing finished at: $(date)"

# cleanup
rm -r TMP_EBB
rm -r TMP_EPICDB
rm -r scores
rm -r EPICDB.BLOOD.METHYL.HERIT
