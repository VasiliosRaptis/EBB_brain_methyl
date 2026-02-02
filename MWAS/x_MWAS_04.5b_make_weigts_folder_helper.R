library(dplyr);library(data.table)

phe <- fread('/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/tensorqtl/phenotypes.int.bed')
profile <- fread('EBB.BRAIN.METHYL.HERIT.profile')
pos <- 
phe[,c(1,2,3,4)] %>% 
  filter(name %in% profile$id) %>% 
  dplyr::rename(CHR = "#chr", PO = start, P1 = end, ID = name) %>% 
  mutate(CHR = gsub('chr', "", CHR))
 
wgt <- 
fread('EBB.BRAIN.METHYL.HERIT.list', header =F) %>% 
  dplyr::rename(WGT=V1) %>% 
  mutate(
    ID = gsub('EBB.BRAIN.METHYL.HERIT/', "", WGT), 
    ID = gsub('.wgt.RDat', "", ID)
    )

pos2 <- inner_join(wgt, pos, by='ID')
fwrite(pos2, file='EBB.BRAIN.METHYL.HERIT.pos', col.names=T, row.names=F, quote=F, sep="\t")
