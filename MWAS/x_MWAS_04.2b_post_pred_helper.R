# --- 
# get the heritable CpG sites and create a list with the paths to their corresponding weights  , to be used by the FUSION.profile_wgt.R script
# Inputs: (1) the .profile file for all CpGs (2) the name of the list with the heritable CpGs weights paths (output name)
# ---

# libraries
library(data.table)
library(dplyr)
# wd 
setwd('/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/MWAS/')

# input / output
arg = commandArgs(trailingOnly=T)
input  = arg[1]
output = arg[2]

# load wgt profile (all CpGs)
profile.all <- fread(input)
profile.hsq <- profile.all %>% filter(hsq.pv < 0.05 & hsq < 1 & hsq > 0)

# compare hsq with previous GCTA estimation
#a <- fread('../heritability/results/h2results.txt')
#test <- inner_join(profile.hsq, a, by = c('id'='cpg'))
#cor(test$hsq, test$h2)
#plot(test$hsq, test$h2); abline(a=0,b=1, col='red')

## summary of models in hsq significant CpGs:

# % best model
cols <- c("top1.pv", "enet.pv", "lasso.pv")
#cols <- c("top1.r2", "enet.r2", "lasso.r2")

best_model <-
profile.hsq %>%
  filter(!is.na(enet.pv) & !is.na(lasso.pv)) %>% 
  mutate(
    best_model = cols[
      max.col(-select(., all_of(cols)), ties.method = "first") # lowest pvalue
    ]
  ) %>%
  # filter r2 > 0.01
  mutate(best_r2 = gsub('.pv', '.r2', best_model)) %>%
  rowwise() %>%
  mutate(best_model_r2 = cur_data()[[best_r2]]) %>%
  ungroup() %>%
  filter(best_model_r2 > 0.01)

print('% model is best (when hsq p < .05):')
best_model %>%  
  group_by(best_model) %>%
  summarise(pct = 100*(n()/nrow(.)))

# avarage r2 
print('avarage r2 (when hsq p < .05):')
best_model %>% 
  summarise(across(c(top1.r2, enet.r2, lasso.r2), ~ round(mean(.x, na.rm = TRUE), 4)))

# make list with hsq significant CpGs - output
list.hsq <-
best_model %>% 
  mutate(list = paste0('/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/MWAS/OUT/', id, '.wgt.RDat')) %>% 
  select(list)
  
print(paste('# of CpG weights with hsq p < .05 and r2 > 0.01:', nrow(list.hsq)))

write.table(list.hsq, output, col.names=F, row.names=F, quote=F)

