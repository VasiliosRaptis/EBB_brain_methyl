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
#expression_bed = rootdir + 'tensorqtl/phenotypes.bed'
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

### Preprocessing

# make dummy variable for categorical covariates
covariates_df = covariates_df.assign(tissue = [1 if tissue == 'Cortex' else 0 for tissue in covariates_df['tissue']])
covariates_df = covariates_df.assign(sex = [1 if sex == 'M' else 0 for sex in covariates_df['sex']])
covariates_df = covariates_df.assign(batch = [1 if batch == '2056G' else 0 for batch in covariates_df['batch']])

# fix chrom name in variant_df: 'chr[chr]'
variant_df = variant_df.assign(chrom = ['chr' + chrom for chrom in variant_df['chrom']])

### cis-QTL: nominal p-values for all variant-phenotype pairs
window = 500000

# map all cis-associations (results for each chromosome are written to file)

# all genes
cis.map_nominal(genotype_df, variant_df,
                phenotype_df,
                phenotype_pos_df,
                prefix, 
                covariates_df=covariates_df,
                window=window,
                output_dir=outdir)
# genes on chr22 ~ 1min
#cis.map_nominal(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                prefix, 
#                covariates_df=covariates_df,
#                window=window,
#                output_dir=outdir)

# load results
#pairs_df = pd.read_parquet(f'{outdir}/{prefix}.cis_qtl_pairs.chr22.parquet')
#pairs_df.head()

### cis-QTL: empirical p-values for phenotypes
window = 500000
nperm  = 10000

# all genes
cis_df = cis.map_cis(genotype_df, variant_df,
                phenotype_df,
                phenotype_pos_df,
                covariates_df=covariates_df,
                window=window,
                nperm=nperm,
                seed=123456)

# genes on chr22 ~ 18.6 min
#cis_df = cis.map_cis(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                covariates_df=covariates_df,
#                window=window,
#                nperm=nperm,
#                seed=123456)

# compute q-values (in practice, this must be run on all genes, not a subset)
#cis_df = pd.read_parquet(f'{outdir}/{prefix}_cis_perm.parquet')
#post.calculate_qvalues(cis_df, fdr=0.05, qvalue_lambda=0.85)

# save
cis_df_out = outdir + "/" + prefix + ".cis_qtl_perm_all.parquet"
cis_df_out_csv = outdir + "/" + prefix + ".cis_qtl_perm_all.csv"

cis_df.to_parquet(cis_df_out)
cis_df.to_csv(cis_df_out_csv)

### calcuate q-values (run the R script: tensorqtl_02_qvalues.R)


### test shuffled ids
#from random import shuffle
#
#rand_samples = list(phenotype_df.columns.values)
#shuffle(rand_samples)
#
#phenotype_df.columns = rand_samples
#covariates_df.index = rand_samples
#
#cis_df_rand22 = cis.map_cis(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                covariates_df=covariates_df,
#                window=window,
#                nperm=nperm,
#                seed=123456)
#
#cis_df_rand22.to_csv('test_perm.csv')
#
#cis.map_nominal(genotype_df, variant_df,
#                phenotype_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                phenotype_pos_df.loc[phenotype_pos_df['chr'] == 'chr22'],
#                'test', 
#                covariates_df=covariates_df,
#                window=window,
#                output_dir=outdir)
#



