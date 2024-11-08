"""
Combines selected regulons with all their previously computed metrics into a single dataframe.
"""

import pandas as pd
import sys

CLUSTER_ID_FILE         = sys.argv[1]  
REGULONS_FILE           = sys.argv[2]
ALIASE_FILE             = sys.argv[3]
REGULON_ENRICHMENT_FILE = sys.argv[4]
CENTRALITIES_FILE       = sys.argv[5]
ALL_MARKERS_FILE        = sys.argv[6]
OUTPUT_FILE             = sys.argv[7]


# ======== dataframe that will collect all the metrics computed for the selected regulons ========
ranking_df = pd.DataFrame()


# ======== load regulons filtered according to parameters specified by the user ========
regulons_df = pd.read_csv(REGULONS_FILE, sep='\t', names=['TF', 'Cluster', 'TGs'])  # here, clusters are identified as "Cluster_X"
regulons_df['Regulon'] = regulons_df['TF'] + '_' + regulons_df['Cluster']  # combine TF and cluster names to create a unique regulon name
regulons_df['#TGs'] = regulons_df['TGs'].str.split(',').str.len()  # count TGs of each regulon
regulons_dict = regulons_df.set_index('Regulon')['#TGs'].to_dict()
# collect the relevant data
ranking_df = regulons_df[['TF', 'Cluster', '#TGs']].copy()  # collect regulon identifiers + #TGs per regulon
regulons_per_cluster_dict = regulons_df.groupby('Cluster')['Regulon'].count().to_dict()
ranking_df['totRegInCluster'] = ranking_df['Cluster'].apply(lambda x: regulons_per_cluster_dict[x])


# ======== load gene aliases for the retrieved TFs ========
tf_alias_df = pd.read_csv(ALIASE_FILE, sep='\t') # columns: 'locus_name' - 'symbol' - 'full_name'
tf_alias_df.columns = ['TF', 'alias', 'description']
tf_alias_df = tf_alias_df.fillna("")
# add entries for the missing TFs, using TF identifiers for the three columns
missing_tfs = set(ranking_df['TF'].unique()) - set(tf_alias_df['TF'])
tf_alias_df = pd.concat([tf_alias_df, pd.DataFrame({ 'TF': list(missing_tfs), 'alias': list(missing_tfs), 'description': list(missing_tfs)})], ignore_index=True)
tf_alias_dict = tf_alias_df.groupby('TF')['alias'].agg('/ '.join).to_dict()
# collect the relevant data
ranking_df['alias'] = ranking_df['TF'].apply(lambda x: tf_alias_dict[x])
ranking_df['clusterId'] = ranking_df['Cluster'].apply(lambda x: x.split('_')[1])


# ======== load cluster annotations ========
cluster_id_annotation_df = pd.read_csv(CLUSTER_ID_FILE, sep='\t', names=['clusterId', 'celltype'], dtype={'clusterId': 'str'})
# collect the relevant data
ranking_df = ranking_df.merge(cluster_id_annotation_df, on='clusterId', how='left')  # add 'celltype' column
ranking_df['cluster'] = ranking_df['celltype'] + "_" + ranking_df['Cluster']


# ======== load cluster specificity for the list of regulons retrieved previously ========
regulon_enrichment_df = pd.read_csv(REGULON_ENRICHMENT_FILE, sep='\t', comment='#')
regulon_enrichment_df = regulon_enrichment_df.rename(columns={'set_id': 'TF', 'ftr_id': 'Cluster', 'q-val': 'qval_cluster'})
regulon_enrichment_df['Regulon'] = regulon_enrichment_df['TF'] + '_' + regulon_enrichment_df['Cluster']  # combine TF and cluster names to create a unique regulon name
regulon_enrichment_df = regulon_enrichment_df[['TF', 'Cluster', 'qval_cluster']]
# collect the relevant data
ranking_df = ranking_df.merge(regulon_enrichment_df, on=['TF', 'Cluster'], how='left')  # add 'qval_cluster' column


# ======== load network centrality measures ========
# the original dataframe contains TFs as rows, and, for each cluster, three columns: 'degout_Cluster_X', 'clos_Cluster_X' and 'bet_Cluster_X'
centrality_orig_df = pd.read_csv(CENTRALITIES_FILE, sep='\t', index_col=0)
centrality_df = regulons_df[['TF', 'Cluster']].copy()  # retrieve identifiers of the selected regulons
centrality_df['out-degree'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"degout_{row['Cluster']}"], axis=1)
centrality_df['closeness'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"clos_{row['Cluster']}"], axis=1)
centrality_df['betweenness'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"bet_{row['Cluster']}"], axis=1)
# collect the relevant data
ranking_df = ranking_df.merge(centrality_df, on=['TF', 'Cluster'], how='left')  # add 'out-degree', 'closeness' and 'betweenness' columns


# ======== load cluster DEGs ========
# this is needed to check whether TFs are DE or not
markers_df = pd.read_csv(ALL_MARKERS_FILE, sep='\t',dtype={'cluster': str})
markers_df = markers_df.rename(columns={'gene': 'TF', 'cluster': 'clusterId'})
markers_df = markers_df[markers_df['p_val_adj'] <= 0.05]
# collect the relevant data
ranking_df['isTF_DE'] = ranking_df.merge(markers_df, on=['TF', 'clusterId'], how='left', indicator=True)['_merge'].map({'both': 1, 'left_only': 0}).astype(int)


# ======== add GO-related information ========
# as this is a standard procedure, GO annotations aren't used --> add default values
ranking_df['hasTFrelevantGOterm'] = "unknown_TF"
ranking_df['GOterm'] = "-"
ranking_df['GOdescription'] = "-"


# ======== select final columns and reorder them ========
ranking_df = ranking_df[['TF', 'alias', 'hasTFrelevantGOterm', 'GOterm', 'GOdescription', 'cluster', 'celltype', 'isTF_DE',
                         'totRegInCluster', '#TGs', 'qval_cluster', 'out-degree', 'closeness', 'betweenness']]


# ======== save the resulting dataframe ========
ranking_df.to_csv(OUTPUT_FILE, sep='\t', index=None)