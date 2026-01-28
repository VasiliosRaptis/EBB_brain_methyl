#!/bin/bash

#$ -cwd
#$ -N wgt
#$ -o ./qsub_logs/wgt3.$JOB_ID.$TASK_ID.o
#$ -e ./qsub_logs/wgt3.$JOB_ID.$TASK_ID.e
#$ -l h_rt=47:59:00
#$ -l h_rss=5G
##$ -pe sharedmem 4
#$ -t 14-24
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
cd /exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/MWAS

module load igmm/apps/plink/1.90b7.2
module load igmm/apps/plink/2.00a6LM

software=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/software
gcta_nr_robust=${software}/fusion_twas-master/gcta_nr_robust
gcta=${software}/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64 # newer version
gemma=${software}/gemma/gemma-0.98.5-linux-static-AMD64
ldref=${software}/fusion_twas-master/LDREF # in hg37 positions

EBB=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation
WD=$EBB/MWAS
TMP=$WD/TMP
INP=$WD/INP
OUT=$WD/OUT
#OUT=$WD/OUT_BSLMM


mkdir -p $TMP
mkdir -p $INP
mkdir -p $OUT

### extract SNPs in LDREF - RUN ONCE
# pgen_prefix=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata/plink_files/pgen/imputed_allchr_newIDs_noOutlier

# cat $ldref/1000G.EUR.*.bim  | awk '{print $1":"$4":"$5":"$6}' > $INP/ldref.extract1 # ref:alt
# cat $ldref/1000G.EUR.*.bim  | awk '{print $1":"$4":"$6":"$5}' > $INP/ldref.extract2 # alr:ref
# cat $EBB/tensorqtl/results/EBB.cis_qtl_indep.csv | awk -F ',' 'NR>1{print $8}' > $INP/meqtls.extract3
# cat $INP/ldref.extract1 $INP/ldref.extract2 $INP/meqtls.extract3 > $INP/ldref.extract

# plink2 --pfile $pgen_prefix --extract $INP/ldref.extract --make-bed --out $INP/ldref

### split phenotype file into chunks - RUN ONCE
# split -n l/50 --numeric-suffixes=1 $EBB/tensorqtl/phenotypes.int.bed $TMP/temp_pheno_bed_
# sed -i '1d' $TMP/temp_pheno_bed_01 # remove header

### get sample names - RUN ONCE
# head -n1 $EBB/tensorqtl/phenotypes.int.bed | tr '\t' '\n' | tail -n+5 | awk '{ print $1,$1 }' > $INP/samples.id

### covariate file - RUN ONCE
#awk '{print $1,$2,$3,$4}' $EBB/heritability/input/qcovar_for_gcta.txt > $INP/covar.txt # sex, age only
#head -n1 $EBB/tensorqtl/covariates.txt | awk '{print $1,$1,$2,$3,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23}' | cat - $EBB/heritability/input/qcovar_for_gcta.txt > $INP/covar.txt # all covariates

### BATCH STARTS HERE

pheno_bed_chunk=$TMP/temp_pheno_bed_$chunk

a=($(wc -l $pheno_bed_chunk))
n_cpgs=${a[0]} # total no. of cpgs 

echo "Analysis for chunk $chunk started at: $(date)"

### loop through each cpg of chunk file
for i in $(seq 1 "$n_cpgs"); do
  
  start=$(date +%s)
  # get cpg position
  cpg_name=$(awk -v i=$i 'NR==i {print $4}' $pheno_bed_chunk)
  cpg_chr=$(awk -v i=$i 'NR==i {print $1}' $pheno_bed_chunk | sed 's/chr//g')
  cpg_pos=$(awk -v i=$i 'NR==i {print $3}' $pheno_bed_chunk)
  
  # get cis-window
  wind=500000
  cpg_pos1=$((cpg_pos + wind)) 
  cpg_pos0=$((cpg_pos - wind)) 
  if (( cpg_pos0 < 0 )); then
      cpg_pos0=0
  fi
  
  # echo $cpg_name $cpg_chr $cpg_pos0 $cpg_pos1
  
  # pull phenotype for cpg i and make pheno file for plink
  awk -v i=$i 'NR==i' $pheno_bed_chunk | tr '\t' '\n' | tail -n+5 | paste -d" " $INP/samples.id - > $TMP/$cpg_name.pheno
  
  # Get the cis locus genotypes for all samples and set current gene expression as the phenotype
  plink --bfile $INP/ldref --chr $cpg_chr --from-bp $cpg_pos0 --to-bp $cpg_pos1 --pheno $TMP/$cpg_name.pheno --make-bed --out $TMP/$cpg_name --threads $(nproc) &>/dev/null 
  
  # calculate prediction weights
  # IMPORTANT: for models without bslmm use: --out $OUT/$cpg_name
  # IMPORTANT: for models with bslmm use: --out $cpg_name  & TMP/.* (not a variable)
  # covariates: phenotype ~ covariate residuals are scaled 
  Rscript ${software}/fusion_twas-master/FUSION.compute_weights.R \
    --bfile $TMP/$cpg_name \
    --covar $INP/covar.txt \
    --tmp $TMP/$chunk_$cpg_name.tmp \
    --out $OUT/$cpg_name \
    --PATH_gcta $gcta \
    --PATH_gemma $gemma \
    --models lasso,top1,enet \
    --hsq_p 1 &>/dev/null 

  #mv $cpg_name.wgt.RDat $OUT/
  
  # cleanup
  rm $TMP/$cpg_name*
  
  end=$(date +%s)
  echo "finished cpg $cpg_name, $i out of $n_cpgs in $((end - start)) seconds"
done

echo "Analysis for chunk $chunk finised at: $(date)"


