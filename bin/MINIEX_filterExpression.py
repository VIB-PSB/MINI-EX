"""
Filters the regulons based on two criteria:
- expression percentage (in each regulon, its TF must be expressed in at least X percent of cells of its cluster)
- target gene cluster specificity (in each regulon, only target genes that are DE in its cluster are kept)
"""

import collections
import pandas as pd
import sys

EXPRESSION_MATRIX_FILE  = sys.argv[1]
TF_LIST_FILE            = sys.argv[2]
CELLS_TO_CLUSTERS_FILE  = sys.argv[3]
EXPRESSION_PERCENTAGE   = sys.argv[4]
CLUSTER_ENRICHMENT_FILE = sys.argv[5]
OUTPUT_FILE             = sys.argv[6]


cells_to_clusters_df = pd.read_csv(CELLS_TO_CLUSTERS_FILE, sep='\t', header=None, names=['cell', 'cluster'], dtype={'cluster': str})
cells_to_clusters_dict = cells_to_clusters_df.set_index('cell')['cluster'].to_dict()
# how many cells in each cluster; ex: Counter({'clusterA': 150, 'clusterB': 728})
cells_to_clusters_counter = collections.Counter(list(cells_to_clusters_dict.values()))

tf_list = pd.read_csv(TF_LIST_FILE, sep='\t', header=None)[0].to_list()

expression_matrix_df = pd.read_csv(EXPRESSION_MATRIX_FILE, sep='\t', index_col=0)
expression_matrix_df = expression_matrix_df[expression_matrix_df.index.isin(tf_list)]  # subset matrix to only keep TFs
expression_matrix_df = expression_matrix_df.transpose()
expression_matrix_df = expression_matrix_df[expression_matrix_df.index.isin(cells_to_clusters_dict.keys())]  # subset matrix to only keep annotated cells
expressed_tfs = list(expression_matrix_df.columns)

# identify valid clusters for each expressed TF (valid = TF is expressed in the specified percentage of cells of that cluster)
tf_to_valid_clusters_dict = {}
for tf in expressed_tfs:
    cells_where_tf_is_expressed = expression_matrix_df[expression_matrix_df[tf] > 1].index.tolist()
    # in how many cells the TF is expressed for each cluster; ex: Counter({'clusterA': 14, 'clusterB': 240})
    cells_per_cluster = collections.Counter(cells_to_clusters_dict[cell] for cell in cells_where_tf_is_expressed)
    # filter clusters based on the expression percentage threshold
    threshold = lambda cluster: cells_per_cluster[cluster] >= round((int(EXPRESSION_PERCENTAGE) / 100) * cells_to_clusters_counter[cluster])
    valid_clusters = [cluster for cluster in cells_per_cluster if threshold(cluster)]
    tf_to_valid_clusters_dict[tf] = valid_clusters if valid_clusters else 'none'

# retrieve target genes for pairs (TF, cluster) identified as valid
cluster_enrichment_df = pd.read_csv(CLUSTER_ENRICHMENT_FILE, sep='\t', comment='#')
cluster_enrichment_df = cluster_enrichment_df.rename(columns={'set_id': 'TF', 'ftr_id': 'cluster', 'hits': 'target_genes'})
cluster_enrichment_df['cluster_id'] = cluster_enrichment_df['cluster'].str.split('_', 1).str[1]
cluster_enrichment_df = cluster_enrichment_df[cluster_enrichment_df.apply(lambda row: row['cluster_id'] in tf_to_valid_clusters_dict.get(row['TF'], []), axis=1)]
cluster_enrichment_df = cluster_enrichment_df[['TF', 'cluster', 'target_genes']]
cluster_enrichment_df.to_csv(OUTPUT_FILE, sep='\t', index=None, header=False)