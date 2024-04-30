import numpy as np
import scanpy as sc
import pandas as pd
from bioservices import BioMart

# load ependymal count data
data_path = '../../data/ependymal/non_neuronal.h5ad'
data = sc.read_h5ad(data_path)

# get homo sapiens pgenes from BioMart
server = BioMart()
xml_query = """
    <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
        <Dataset name="hsapiens_gene_ensembl" interface="default">
            <Attribute name="ensembl_gene_id" />
            <Attribute name="gene_biotype" />
        </Dataset>
    </Query>
"""
response = server.query(xml_query)

# convert to dataframe
protein_coding_genes_df = pd.DataFrame([x.split('\t') for x in response.strip().split('\n')],
                                        columns=['ensembl_gene_id', 'gene_biotype'])

# remove non protein-coding genes
protein_coding_genes_df = protein_coding_genes_df[protein_coding_genes_df['gene_biotype'] == 'protein_coding']
protein_coding_mask = data.var.index.isin(protein_coding_genes_df['ensembl_gene_id'])
data = data[:, protein_coding_mask]

# normalize counts
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)

# rank genes by mean expression
data.var['mean_normalized_expression'] = np.mean(data.X.toarray(), axis=0)
genes = data.var.copy(deep=True)
genes.sort_values(by='mean_normalized_expression', ascending=False, inplace=True)

# cleanup
genes['Rank'] = np.arange(1, len(genes) + 1)
genes.set_index('Gene', inplace=True)
genes.drop(columns=['Biotype', 'End', 'Start', 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype'], inplace=True)

# save mean gene expression
genes.to_csv('./ependymal_mean_expression.csv')