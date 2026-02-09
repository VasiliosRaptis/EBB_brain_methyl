### in R:
setwd("/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/genotypingdata")
library(dplyr)
library(data.table)
library(ggplot2)

## load king matrix and flatten
king <- fread('plink_files/rel_mat/non_imputed_postqc.king') %>% as.matrix
id <- fread('plink_files/rel_mat/non_imputed_postqc.king.id')
colnames(king) <- id$IID
rownames(king) <- id$IID

flattenKingMatrix <- function(cormat) {
  lt <- lower.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[lt]],
    column = rownames(cormat)[col(cormat)[lt]],
    coef  = (cormat)[lt]
    )
}

king_df <- flattenKingMatrix(king)

## summary (no relatives)
summary(king_df$coef)

## plot king coefficients
png('plink_files/rel_mat/king_plot.png')
hist(king_df$coef, xlim=c(-.05,0.5))
dev.off()

## how many PCs to include?
eigenval <- fread('plink_files/rel_mat/non_imputed_postqc_pca.eigenval')

## make scree plot -> eigenvalues flatten at PC=4
png('plink_files/rel_mat/eigenval_plot.png')
plot(eigenval$V1, main='GRM Scree plot')
dev.off()



## find population outliers (+-5df from POP PC1&2 mean)
pca <- fread('plink_files/rel_mat/non_imputed_postqc_pca.eigenvec')

pc1_min <- mean(pca$PC1) - sd(pca$PC1)*5
pc1_max <- mean(pca$PC1) + sd(pca$PC1)*5
pc2_min <- mean(pca$PC2) - sd(pca$PC2)*5
pc2_max <- mean(pca$PC2) + sd(pca$PC2)*5

pca2 <-
pca %>%
  mutate(outlier = ifelse( PC1 < pc1_min | PC1 > pc1_max | PC2 < pc2_min | PC2 > pc2_max, 1, 0)) %>%
  mutate(outlier = as.factor(outlier))

png('plink_files/rel_mat/pca_plot.png')
ggplot(pca2, aes(x=PC1, y=PC2, color=outlier)) + geom_point() + theme_bw()
dev.off()

## "SD037/14B" is population outlier to remove
outlier.ids <- pca2 %>% filter(outlier == 1) %>% pull(IID)

