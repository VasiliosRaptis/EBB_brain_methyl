import pandas as pd
import torch
import tensorqtl
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
outdir = rootdir + 'tensorqtl/results'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)#.T

# PLINK reader for genotypes
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df

# load cis-qtl results
cis_df2   = pd.read_csv('/exports/cmvm/eddie/smgphs/groups/Quantgen/Users/vasilis/PHD/EBB_methylation/tensorqtl/results/EBB.cis_qtl_perm_all_qval.csv')  # output from map_cis + qvalues
cis_df2.set_index('phenotype_id', inplace=True)

### Preprocessing

# make dummy variable for categorical covariates
covariates_df = covariates_df.assign(tissue = [1 if tissue == 'Cortex' else 0 for tissue in covariates_df['tissue']])
covariates_df = covariates_df.assign(sex = [1 if sex == 'M' else 0 for sex in covariates_df['sex']])
covariates_df = covariates_df.assign(batch = [1 if batch == '2056G' else 0 for batch in covariates_df['batch']])

# fix chrom name in variant_df: 'chr[chr]'
variant_df = variant_df.assign(chrom = ['chr' + chrom for chrom in variant_df['chrom']])

### conditionally independent cis-qtl
window=500000
indep_df = cis.map_independent(genotype_df, variant_df, cis_df2, phenotype_df, phenotype_pos_df, covariates_df, window=window)

# save
indep_outdir = outdir + '/EBB.cis_qtl_indep.csv'
indep_df.to_csv(indep_outdir)



