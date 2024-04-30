import numpy as np
import pandas as pd

risk_loci = {}

# Alzheimer's
AD_df = pd.read_excel('./AD_gwas_S5.xlsx', sheet_name='Supplementary Table 5')['Gened']

AD_df.dropna(inplace=True)
risk_loci["AD"] = AD_df.unique().tolist()

# Multiple sclerosis
MS_df = pd.read_excel('./MS_gwas_S18.xlsx', sheet_name='ST18')

MS_genes = []
for col in ['Exonic genes', 'eQTL genes', 'Regulatory network']:
    gs = MS_df[col].dropna().tolist()
    for g in gs:
        MS_genes.extend(g.split('|'))

risk_loci['MS'] = np.unique(MS_genes).tolist()

# Schizophrenia
risk_loci['SCZ'] = pd.read_excel('./SCZ_gwas_S12.xlsx', sheet_name='Prioritised')['Symbol.ID'].tolist()

# Epilepsy
risk_loci['Epilepsy'] = ['FANCL', 'BCL11A', 'SCN3A', 'SCN2A', 'TTC21B',
                        'SCN1A', 'HEATR3', 'BRD7', 'STAT4', 'PCDH7',
                        'GABRA2', 'KCNN2', 'ATXN1', 'PNPO', 'GRIK1',
                        'STX1B', 'ZEB2', 'C3ORF33', 'SLC33A1', 'KCNAB1', 'GJA1']

# Stroke
Stroke_df = pd.read_excel('./Stroke_gwas_S4.xlsx', sheet_name='Suppl.Table16', skiprows=4, usecols='A:K')

risk_loci['Stroke'] = Stroke_df['Gene'].dropna().tolist()

# Hydroencephalus
Hydro_df = pd.read_excel('./Hydro_gwas.xlsx', header=None)

risk_loci['Hydroencephalus'] = Hydro_df[0].dropna().tolist()

# combine results
ependymal_gex = pd.read_csv('./ependymal_median_expression.csv')

for dx, genes in risk_loci.items():
    ependymal_gex[dx] = ependymal_gex['Gene'].isin(genes)

# remove genes with no overlap
nonzero_mask = ependymal_gex.loc[:,'AD':'Stroke'].sum(axis=1) > 0
ependymal_gex_masked = ependymal_gex[nonzero_mask]
ependymal_gex_masked['count'] = ependymal_gex_masked.loc[:,'AD':'Hydroencephalus'].sum(axis=1).values

ependymal_gex_masked.to_csv('./ependymal_gex_gwas_overlap.csv', index=False)


