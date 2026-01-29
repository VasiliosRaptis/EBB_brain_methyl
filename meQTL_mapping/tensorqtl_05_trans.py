import pandas as pd
import torch
import tensorqtl
import sys
from tensorqtl import pgen, cis, trans, post
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")

### Setup

# define paths to data
rootdir = '/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/'
plink_prefix_path =  rootdir + 'genotypingdata/plink_files/pgen/imputed_allchr_newIDs_noOutlier'
expression_bed = rootdir + 'tensorqtl/phenotypes.int.bed'
covariates_file = rootdir + 'tensorqtl/covariates.txt'
prefix = 'EBB'
outdir = rootdir + 'tensorqtl/results/trans'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)#.T

# PLINK reader for genotypes
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df

### Preprocessing

# make dummy variable for categorical covariates
covariates_df = covariates_df.assign(tissue = [1 if tissue == 'Cortex' else 0 for tissue in covariates_df['tissue']])
covariates_df = covariates_df.assign(sex = [1 if sex == 'M' else 0 for sex in covariates_df['sex']])
covariates_df = covariates_df.assign(batch = [1 if batch == '2056G' else 0 for batch in covariates_df['batch']])

# fix chrom name in variant_df: 'chr[chr]'
variant_df = variant_df.assign(chrom = ['chr' + chrom for chrom in variant_df['chrom']])

### trans qtl mapping - per chromosome
chrom = 'chr' + str(sys.argv[1]) # pass chrom number from bash script

trans_df = trans.map_trans(genotype_df, phenotype_df.loc[phenotype_pos_df['chr'] == chrom], covariates_df,
                           return_sparse=True, pval_threshold=1e-5, maf_threshold=0.01,
                           batch_size=20000)
# remove cis-associations
window = 500000
trans_df = trans.filter_cis(trans_df, phenotype_pos_df, variant_df, window=window)

# save
trans_outdir = outdir + f"/EBB.trans_qtl_{chrom}.parquet"
trans_df.to_parquet(trans_outdir)

