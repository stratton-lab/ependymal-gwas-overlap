import numpy as np
import pandas as pd
import scanpy as sc

data_path = ''
data = sc.read_h5ad(data_path)
gwas_genes = pd.read_csv('./ependymal_gex_gwas_overlap.csv', index_col=0)

celltypes = data.obs['cell_type'].unique().tolist()

# normalize count data
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)

n_top = 100
mean_gex_df = pd.DataFrame(np.zeros((gwas_genes.shape[0], len(celltypes))), index=gwas_genes.index, columns=celltypes)

# calculate mean expression across all celltypes
for i in range(n_top):
    print(i)
    mean_gex = []

    goi = gwas_genes.index[i]
    for ct in celltypes:
        mean_gex.append(np.mean(data.X[data.obs['cell_type'] == ct,data.var['Gene']==goi]))

    mean_gex_df.iloc[i,:] = mean_gex

df_merged = pd.merge(gwas_genes, mean_gex_df, left_index=True, right_index=True)

df_merged.to_csv('./ependymal_gwas_genes_compared.csv')