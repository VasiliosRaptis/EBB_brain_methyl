# --- 
# computes the actual CpG level ~ imputed CpG level correlation 
# Inputs: (1) the .sscore plink file for a CpG phenotype and (2) the CpG name
# ---

arg = commandArgs(trailingOnly=T)

# load score file
scorefile = arg[1]
score <- read.table(scorefile, header=T, comment.char="")

# get model used
mod <- gsub("_AVG", "", names(score)[6])

# get CpG name
cpgname = arg[2]

# pearson's correlation 
cpg.cor.test<- cor.test(score$PHENO1, score[,6], method='p')
cpg.cor     <- cpg.cor.test$estimate
cpg.cor.pv  <- cpg.cor.test$p.value
cpg.cor.lci <- cpg.cor.test$conf.int[1]
cpg.cor.uci <- cpg.cor.test$conf.int[2]

r2 <- cpg.cor^2

# write output
out = cbind(cpgname, mod, cpg.cor, cpg.cor.pv, cpg.cor.lci, cpg.cor.uci, r2)

write.table(out, quote=F , row.names=F , col.names=F , sep='\t' )
