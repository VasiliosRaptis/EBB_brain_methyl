#!/bin/bash

#$ -cwd
#$ -N int
#$ -o ./qsub_logs/int.o
#$ -e ./qsub_logs/int.e
#$ -l h_rt=23:00:00
#$ -l h_vmem=8G
#$ -pe sharedmem 8
##$ -t 1-22
##$ -tc 22

## initialise the modules framework
source /etc/profile.d/modules.sh
#module load python/3.12.9
module load anaconda/2024.02

## run int script
conda activate jpt
Rscript meffil_00_int.R