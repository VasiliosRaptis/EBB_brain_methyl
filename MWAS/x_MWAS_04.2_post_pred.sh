#!/bin/bash

#$ -cwd
#$ -N wgtlist
#$ -o ./qsub_logs/wgtlist.$JOB_ID.$TASK_ID.o
#$ -e ./qsub_logs/wgtlist.$JOB_ID.$TASK_ID.e
#$ -l h_rt=35:00:00
#$ -l h_rss=4G
#$ -pe sharedmem 4
##$ -t 1-50
##$ -tc 50

source /etc/profile.d/modules.sh

software=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/software
EBB=/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation
name=EBB.BRAIN.METHYL


### make weights summaries for ALL CpGs
cd /exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/MWAS/OUT
# get list of all CpGs (list of paths)
find . -name '*wgt*' -exec realpath {} + > ../$name.ALL.list
# make profile table & summaries
cd /exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/MWAS/
Rscript $software/fusion_twas-master/utils/FUSION.profile_wgt.R $name.ALL.list > $name.ALL.profile 2> $name.ALL.profile.err

### make weights summaries for HERITABLE CpGs (hsq p < 0.05 & r2 > 0.01)
cd /exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/MWAS
# get list of heritable CpGs (list of paths)
Rscript $EBB/scripts/x_MWAS_04.2b_post_pred_helper.R $name.ALL.profile $name.HERIT.list
# make profile table & summaries
Rscript $software/fusion_twas-master/utils/FUSION.profile_wgt.R $name.HERIT.list > $name.HERIT.profile 2> $name.HERIT.profile.err
