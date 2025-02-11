"""
Computes network centrality metrics (out degree, closeness and betweenness) for each regulon.
In addition, closeness and betweenness are also computed for each TF from the original graph
(after the motif filtering step) and stored in the columns "orig_closeness" and "orig_betweenness".
The resulting dataframe contains one row per TF and, for each TF,:
- three columns per cluster (one for each per cluster metric)
- columns "orig_closeness" and "orig_betweenness".
EX: for a dataset with two clusters:
TF    out-degree_Cl0  closeness_Cl0  betweenness_Cl0  out-degree_Cl1  closeness_Cl1  betweenness_Cl1  orig_closeness  orig_betweenness
TF1              0.5            0.2             0.14             0.3            0.4             0.15            0.67               0.8       
TF2              0.2            0.2              0.6            0.58            0.3             0.28            0.43               0.6
"""

import sys
import numpy as np
import pandas as pd
import networkx as nx
import graph_tool as gt

from graph_tool.centrality import closeness, betweenness

       
FINAL_REGULONS_FILE    = sys.argv[1]  # per cluster regulons
ORIGINAL_REGULONS_FILE = sys.argv[2]  # regulons after the motif filtering step
OUPUT_FILE             = sys.argv[3]


# ======== read the per cluster regulons dataframe ========
regulons_df = pd.read_csv(FINAL_REGULONS_FILE, sep='\t', header=None, names=['TF', 'cluster', 'TGs'])
regulons_df['TGs'] = regulons_df['TGs'].str.split(',')  # convert to list of TGs

all_tfs = regulons_df['TF'].unique()
all_clusters = regulons_df['cluster'].unique()


# ======== create dictionary associating clusters to a dictionary associating regulators of that cluster to their TGs ========
# example: { 'Cluster0' : { 'TF1' : ['TG1', 'TG2'], 'TF2' : ['TG2', 'TG3', 'TG4'] } }
cluster_to_tf_to_tg_dict = {
    cluster: regulons_df[regulons_df['cluster'] == cluster]  # filter on each cluster
    .groupby('TF')['TGs']                                    # for each TF of that cluster: collect its target genes
    .apply(lambda tg_lists: list(set(sum(tg_lists, []))))    # flatten and deduplicate
    .to_dict()
    for cluster in all_clusters
}


# ======== prepare the centralities dataframe ========
# the dataframe will contain one row per TF and three columns per cluster (one for each per cluster metric)
# one column for original closeness and one for original betweenness
centrality_metrics = ['out-degree', 'closeness', 'betweenness']
columns = [f"{metric}_{cluster}" for cluster in sorted(all_clusters) for metric in centrality_metrics]  # per cluster metrics
columns += ['orig_closeness', 'orig_betweenness']  # original metrics
centralities_df = pd.DataFrame(float('nan'), index=all_tfs, columns=columns)


# ======== define closeness and betweenness functions ========
def closeness_centrality(graph, number_to_gene):
    
    # create a reversed graph view for in-closeness (to match NetworkX's directionality)
    graph_reversed = gt.GraphView(graph, reversed=True)

    # compute both normalized and raw closeness values to get the number of reachable nodes (n_r)
    closeness_norm = closeness(graph_reversed, norm=True)  # (n_r - 1)/sum_d
    closeness_raw = closeness(graph_reversed, norm=False)  # 1/sum_d

    # calculate n_r - 1 (number of reachable nodes minus 1) for each node
    N = graph_reversed.num_vertices()
    n_r_minus_1 = np.zeros_like(closeness_norm.a)
    valid = ~np.isinf(closeness_raw.a)  # Filter out nodes where closeness_raw is infinity
    n_r_minus_1[valid] = closeness_norm.a[valid] / closeness_raw.a[valid]

    # apply NetworkX's scaling: multiply by (n_r - 1)/(N - 1)
    closeness_scaled = closeness_norm.a * (n_r_minus_1 / (N - 1))

    # Handle isolated nodes (sum_d=0) and NaN values explicitly
    closeness_scaled[np.isinf(closeness_raw.a) | np.isnan(closeness_scaled)] = 0.0

    # convert to dictionary with gene IDs
    closeness_dict = {number_to_gene[i]: closeness_scaled[i] for i in range(graph_reversed.num_vertices())}

    return closeness_dict


def betweenness_centrality(graph, number_to_gene):

    # compute betweenness using Graph-tool
    betweenness_gt = betweenness(graph)

    # convert to dictionary with gene IDs
    betweenness_dict = {number_to_gene[i]: betweenness_gt[0][i] for i in range(graph.num_vertices())}

    return betweenness_dict

# ======== compute per cluster metrics ========
for cluster, tf_to_tg_dict in cluster_to_tf_to_tg_dict.items():

    # create an edge list for the current cluster
    edges = [(tf, tg) for tf, tg_list in tf_to_tg_dict.items() for tg in tg_list]

    # convert genes to numbers (for graph-tool) and store the conversion information for later
    gene_to_number = {node: i for i, node in enumerate(set([node for edge in edges for node in edge]))}
    number_to_gene = {i: node for i, node in enumerate(set([node for edge in edges for node in edge]))}
    edges_numbers = [[gene_to_number[edge[0]], gene_to_number[edge[1]]] for edge in edges]

    # create a directed NetworkX graph
    graph_nx = nx.DiGraph()
    graph_nx.add_edges_from(edges)

    # create a directed Graph-tool graph
    graph_gt = gt.Graph(directed=True)
    graph_gt.add_edge_list(edges_numbers)

    # calculate centralities
    out_degree_dict = nx.out_degree_centrality(graph_nx)
    closeness_dict = closeness_centrality(graph_gt, number_to_gene)
    betweenness_dict = betweenness_centrality(graph_gt, number_to_gene)

    # populate the centralities dataframe
    for tf in all_tfs:
        centralities_df[f'out-degree_{cluster}'][tf] = out_degree_dict.get(tf)
        centralities_df[f'closeness_{cluster}'][tf] = closeness_dict.get(tf)
        centralities_df[f'betweenness_{cluster}'][tf] = betweenness_dict.get(tf)


# ======== compute the original metrics ========
original_regulons_df = pd.read_csv(ORIGINAL_REGULONS_FILE, sep='\t', header=None, usecols=[0, 1], names=['TF', 'TG']).drop_duplicates()

# convert dataframe to edge list
edges = original_regulons_df.values.tolist()

# convert genes to numbers (for graph-tool) and store the conversion information for later
gene_to_number = {node: i for i, node in enumerate(set([node for edge in edges for node in edge]))}
number_to_gene = {i: node for i, node in enumerate(set([node for edge in edges for node in edge]))}
edges_numbers = [[gene_to_number[edge[0]], gene_to_number[edge[1]]] for edge in edges]

# create a directed Graph-tool graph
original_graph_gt = gt.Graph(directed=True)
original_graph_gt.add_edge_list(edges_numbers)

# calculate centralities
closeness_dict = closeness_centrality(original_graph_gt, number_to_gene)
betweenness_dict = betweenness_centrality(original_graph_gt, number_to_gene)

# populate the centralities dataframe
centrality_data = pd.DataFrame({
    'TF': all_tfs,
    'orig_closeness': [closeness_dict.get(tf) for tf in all_tfs],
    'orig_betweenness': [betweenness_dict.get(tf) for tf in all_tfs],
}).set_index('TF')

# add the two columns to the centralities dataframe
centralities_df.update(centrality_data)


# ======== save the resulting dataframe ========
centralities_df.to_csv(OUPUT_FILE, sep='\t')
