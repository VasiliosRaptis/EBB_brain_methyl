#!/bin/bash

### inputs needed: 
## passed by -iin: batch*_predict.txt
## passed by -iin: batch*_predict_blood.txt
## passed by -iin: icd.phe 
## passed by -iin: cov.phe
## passed by -iin: phecode_map.csv
## passed by -iin: mwas_func.py
## passed by -iin: mwas_07.6_bb_MWAS.py
## passed by -iin: mwas_07.6_bb_MWAS.sh
## passed by -icmd: $1 batch

### setup batch
batch=$1

### setup phecode mapping file
PHE_MAP=phecode_map.csv

### setup inputs
ICD=icd.phe 
COV=cov.phe

### setup output file suffix 
OUT=_batch${batch}.csv
    
### setup python
pip install statsmodels
pip install pandas
pip install scipy
pip install numpy

### set timer
start=$(date +%s.%N)
echo "starting batch $batch at $(date)"

### run association analysis
python3 mwas_07.6_bb_MWAS.py $batch $PHE_MAP $ICD $COV $OUT

### print time
end=$(date +%s.%N)
time_i=$(awk "BEGIN { printf \"%.2f\", $end - $start }")
echo "Finished batch $batch in $time_i seconds"

### cleanup