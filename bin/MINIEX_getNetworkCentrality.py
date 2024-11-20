"""
Computes network centrality metrics (out degree, closeness and betweenness) for each regulon.
The resulting dataframe contains one row per TF and three columns per cluster (one for each metric).
"""

import pandas as pd
import networkx as nx
import sys
       
REGULONS_FILE = sys.argv[1]
OUPUT_FILE    = sys.argv[2]


# ======== read the regulons dataframe ========
regulons_df = pd.read_csv(REGULONS_FILE, sep='\t', header=None, names=['TF', 'cluster', 'TGs'])
regulons_df['TGs'] = regulons_df['TGs'].str.split(',')  # convert to list of TGs

all_tfs = regulons_df['TF'].unique()
all_clusters = regulons_df['cluster'].unique()


# ======== create dictionary associating clusters to a dictionary associating regulators of that cluster to their TGs ========
# example: {"Cluster0" : {"TF1" : [TG1, TG2], "TF2" : [TG2, TG3, TG4]} }
cluster_to_tf_to_tg_dict = {
    cluster: regulons_df[regulons_df['cluster'] == cluster]  # filter on each cluster
    .groupby('TF')['TGs']  # for each TF of that cluster: concatenate its target genes
    .apply(lambda tg_lists: list(set(sum(tg_lists, []))))  # flatten and deduplicate
    .to_dict()
    for cluster in all_clusters
}


# ======== prepare the centralities dataframe ========
# the dataframe will contain one row per TF and three columns per cluster (one for each metric)
centrality_metrics = ['out-degree', 'closeness', 'betweenness']
columns = [f"{metric}_{cluster}" for cluster in sorted(all_clusters) for metric in centrality_metrics]
centralities_df = pd.DataFrame(float('nan'), index=all_tfs, columns=columns)

for cluster, tf_to_tg_dict in cluster_to_tf_to_tg_dict.items():
    # create a directed graph for the current cluster
    graph = nx.DiGraph()
    edges = [(tf, tg) for tf, tg_list in tf_to_tg_dict.items() for tg in tg_list]
    graph.add_edges_from(edges)

    # calculate centralities
    out_degree = nx.out_degree_centrality(graph)
    closeness = nx.closeness_centrality(graph)
    betweenness = nx.betweenness_centrality(graph)

    # populate the centralities dataframe
    for tf in tf_to_tg_dict:
        centralities_df[f'out-degree_{cluster}'][tf] = out_degree.get(tf)
        centralities_df[f'closeness_{cluster}'][tf] = closeness.get(tf)
        centralities_df[f'betweenness_{cluster}'][tf] = betweenness.get(tf)


# ======== save the resulting dataframe ========
centralities_df.to_csv(OUPUT_FILE, sep='\t')
