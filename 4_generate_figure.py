import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# load expert annotated GWAS hits (primary and secondary)
data_path = ''
genes_primary = pd.read_excel(data_path, sheet_name='Primary Ependymal Hits')
genes_primary['level'] = 'primary'
genes_secondary = pd.read_excel(data_path, sheet_name='Secondary Ependymal Hits')
genes_secondary['level'] = 'secondary'
genes = pd.concat([genes_primary, genes_secondary], axis=0)

# reorder columns
rearranged_columns = ['Gene', 'Rank', 'MS', 'AD', 'SCZ', 'Epilepsy', 'Stroke', 
                      'Hydrocephalus', 'count', 'ependymal cell', 'choroid plexus epithelial cell', 
                      'oligodendrocyte', 'astrocyte', 'Bergmann glial cell', 
                      'oligodendrocyte precursor cell', 'fibroblast', 'pericyte', 
                      'vascular associated smooth muscle cell', 'endothelial cell', 
                      'CNS macrophage (microglia)', 'neuron', 'leukocyte', 'level']
genes = genes[rearranged_columns]

# organize dataframe
genes.sort_values(by=['MS','AD','SCZ','Epilepsy','Stroke','Hydrocephalus'], ascending=False, inplace=True)
genes.set_index('Gene', inplace=True)
genes.index = genes.index.to_series().str.split(' ', expand=True)[0].values

genes.loc[:,'ependymal cell':'leukocyte'] = genes.loc[:,'ependymal cell':'leukocyte'].values / genes.loc[:,'ependymal cell'].values.reshape(-1,1)

# generate heatmap
matplotlib.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.family': 'Helvetica Neue'})

f, ax = plt.subplots(figsize=(10, 10))

sns.heatmap(genes.loc[:,'ependymal cell':'leukocyte'], annot=False, ax=ax, cmap='GnBu', linewidths=1, linecolor='w', vmax=1)
plt.tick_params(axis='x', which='major', labelsize=10, labelbottom = False, bottom=False, top = False, labeltop=True, rotation=90)

# plt.savefig('./ependymal_heatmap.svg', dpi=600, bbox_inches='tight')