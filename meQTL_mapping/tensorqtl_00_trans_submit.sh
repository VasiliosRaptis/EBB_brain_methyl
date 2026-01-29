#!/bin/bash

#$ -cwd
#$ -N trans
#$ -o ./qsub_logs/trans.$JOB_ID.$TASK_ID.o
#$ -e ./qsub_logs/trans.$JOB_ID.$TASK_ID.e
#$ -l h_rt=47:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
#$ -t 1-22
#$ -tc 22

## initialise the modules framework
source /etc/profile.d/modules.sh
#module load python/3.12.9
module load anaconda/2024.02

## run tensorqtl script
conda activate tensorqtl
python -V
python ./tensorqtl_05_trans.py $SGE_TASK_ID