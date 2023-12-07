from itertools import combinations
from typing import Dict
from scipy import stats
import sys
import pandas as pd


REGULONS_FILE_NAME=sys.argv[1]
OUTPUT_FILE_NAME=sys.argv[2]
PROCEDURE=sys.argv[3]


def run_borda():
    # the regulons dataframe provides all the information about the regulons
    # collected so far, including the values of the individual metrics
    regulons_df = pd.read_csv(REGULONS_FILE_NAME, sep='\t')

    # list of metrics that can be used for borda ranking + their sorting order ('ascending' = lower values have better rank)
    metrics = {'qval_cluster': 'ascending', 'out-degree': 'descending', 'betweenness': 'descending', 'closeness': 'descending', 'GO_enrich_qval': 'ascending'}
    if 'GO_enrich_qval' not in regulons_df.columns:  # if no GO enrichment information is available: remove this metric
        metrics.pop('GO_enrich_qval')

    # compute individual (global and cluster) rankings for each metric and add them as temporary columns to the regulons dataframe
    regulons_df = add_metric_rankings(regulons_df, metrics)

    # depending on the selected procedure, compute Borda values and rankings and add them to the regulons dataframe
    if PROCEDURE == "std":
        regulons_df = run_borda_std(regulons_df, metrics)
    elif PROCEDURE == "ref":
        regulons_df = run_borda_ref(regulons_df, metrics)

    regulons_df = remove_temporary_columns(regulons_df) # remove metric ranks & Borda values => only keep Borda ranks
    regulons_df = regulons_df.sort_values(by="borda_rank")  # sort rows by global Borda ranking

    regulons_df.to_excel(OUTPUT_FILE_NAME, index=None)


def add_metric_rankings(regulons_df: pd.DataFrame, metrics: Dict):
    """
    For each metric, adds two columns to the original dataframe:
    - metricName_metricRank: ranks regulons by sorting the entire dataset
      based on the metric values
    - metricName_metricClusterRank: ranks regulons within individual clusters
      based on the metric values
    Ranking order is defined by the provided metrics dictionary.
    In case of tied values, the highest rank is assigned to all corresponding regulons.

    Parameters
    ----------
    regulons_df : pd.DataFrame
        regulons dataframe
    metrics : Dict
        dictionary of metrics to use, together with their sort order (ascending or descending)

    Returns
    -------
    pd.DataFrame
        original regulon dataframe with additional columns: global and cluster ranking for each metric
    """
    for metric_name, metric_sort_order in metrics.items():
        regulons_df[f'{metric_name}_metricRank'] = regulons_df[metric_name].rank(
            method ='max',
            ascending=True if metric_sort_order == 'ascending' else False,
            na_option='bottom')
        regulons_df[f'{metric_name}_metricClusterRank'] = regulons_df.groupby('cluster')[metric_name].rank(
            method ='max',
            ascending=True if metric_sort_order == 'ascending' else False,
            na_option='bottom')
    return regulons_df


def remove_temporary_columns(regulons_df: pd.DataFrame):
    # removes columns with metric ranks, as well as the computed Borda values and returns the dataframe
    regulons_df = regulons_df.drop(columns=['borda', 'bordaCluster'])
    regulons_df = regulons_df.drop(columns=regulons_df.filter(regex='_metricRank|_metricClusterRank').columns)
    return regulons_df


def run_borda_std(regulons_df: pd.DataFrame, metrics: Dict):
    """
    Computes Borda ranks using the standard approach:
    - for each regulon, its Borda value corresponds to the geometric mean of its metrics
    - the computed Borda values are then ranked in ascending order
    This process is conducted at the global level using global metric ranks
    (resulting in the 'borda_rank' column) and at the per-cluster level
    using per-cluster metric ranks (resulting in the 'borda_clusterRank' column).

    Parameters
    ----------
    regulons_df : pd.DataFrame
        regulons dataframe
    metrics : Dict
        dictionary of metrics to use, together with their sort order (ascending or descending)

    Returns
    -------
    pd.DataFrame
        original dataframe with Borda-related columns
    """
    print("Borda procedure: STANDARD")
    print(f"Default metrics: {', '.join(metrics.keys())}")

    rank_columns = regulons_df.columns[regulons_df.columns.str.endswith('_metricRank')]
    regulons_df['borda'] = stats.gmean(regulons_df[rank_columns], axis=1)
    regulons_df['borda_rank']=regulons_df['borda'].rank(method ='max',ascending=True)

    cluster_rank_columns = regulons_df.columns[regulons_df.columns.str.endswith('_metricClusterRank')]
    regulons_df['bordaCluster'] = stats.gmean(regulons_df[cluster_rank_columns], axis=1)
    regulons_df['borda_clusterRank'] = regulons_df.groupby('cluster')['bordaCluster'].rank(method ='max', ascending=True)
    
    return regulons_df


def run_borda_ref(regulons_df: pd.DataFrame, metrics: Dict):
    """
    Computes Borda ranks using the reference approach:
    - for each metric, compute its weighted ranks by dividing the original ranks by the metric's R50
    - for each combination of metrics (going from a single metric to all the metrics),
      compute its Borda R50, with Borda values defined as a sum of all the weighted
      ranks of the metrics in that combination
    - select the best combination of metrics, defined as having the lowest Borda R50
    - using the selected combination of metrics, compute the global Borda rank (resulting
      in the 'borda_rank' column) and per-cluster Borda rank (resulting in the 'borda_clusterRank'
      column)

    Parameters
    ----------
    regulons_df : pd.DataFrame
        regulons dataframe
    metrics : Dict
        dictionary of metrics to use, together with their sort order (ascending or descending)

    Returns
    -------
    pd.DataFrame
        original dataframe with Borda-related columns
    """
    print("Borda procedure: REFERENCE")

    # compute weighted metric ranks: the initial ranks are divided by the R50 of that metric
    for metric in metrics.keys():
        r50 = get_r50_for_metric(regulons_df, metric)
        regulons_df[f'{metric}_metricRankWeighted'] = regulons_df[f'{metric}_metricRank']/r50
        regulons_df[f'{metric}_metricClusterRankWeighted'] = regulons_df[f'{metric}_metricClusterRank']/r50

    # generate all the combinations of the metrics
    metric_combinations = [comb for i in range(1, len(metrics.keys()) + 1) for comb in combinations(metrics.keys(), i)]

    # compute Borda R50 for each combination of metrics and store it into a dictionary
    metric_combination_r50 = {}
    for metric_combination in metric_combinations:
        weighted_metric_combination = [item + '_metricRankWeighted' for item in metric_combination]
        regulons_df['borda'] = regulons_df[weighted_metric_combination].sum(axis=1)
        regulons_df['borda_metricRank']=regulons_df['borda'].rank(method='max', ascending=True)
        borda_r50 = get_r50_for_metric(regulons_df, 'borda')
        metric_combination_r50[','.join(metric_combination)] = borda_r50

    # identify the best combination of metrics (=having the lowest Borda R50)
    best_metric_combination = min(metric_combination_r50, key=lambda k: metric_combination_r50[k]).split(',')
    print(f"Best metrics: {', '.join(best_metric_combination)}")

    # compute final Borda rank based on the best combination of metrics
    weighted_metric_combination = [item + '_metricRankWeighted' for item in best_metric_combination]
    regulons_df['borda'] = regulons_df[weighted_metric_combination].sum(axis=1)
    regulons_df['borda_rank']=regulons_df['borda'].rank(method='max', ascending=True)
    
    # compute final Borda rank per cluster based on the best combination of metrics
    weighted_cluster_metric_combination = [item + '_metricClusterRankWeighted' for item in best_metric_combination]
    regulons_df['bordaCluster'] = regulons_df[weighted_cluster_metric_combination].sum(axis=1)
    regulons_df['borda_clusterRank'] = regulons_df.groupby('cluster')['bordaCluster'].rank(method ='max', ascending=True)

    return regulons_df


def get_r50_for_metric(regulons_df: pd.DataFrame, metric_name: str):
    """
    Computes the R50 value for the provided metric.
    R50 corresponds to the rank at which half of the TFs
    annotated as relevant are recovered using the provided metric.

    Parameters
    ----------
    regulons_df : pd.DataFrame
        regulons dataframe
    metric_name : str
        name of the metric to use for the computation of R50

    Returns
    -------
    int
        R50 for the provided metric
    """
    sorted_regulons_df=regulons_df.sort_values(by=f"{metric_name}_metricRank")
    relevant_tfs_df = sorted_regulons_df[sorted_regulons_df['hasTFrelevantGOterm'] == 'relevant_known_TF']
    relevant_tfs_to_find = round(len(relevant_tfs_df)/2)
    r50 = relevant_tfs_df[f"{metric_name}_metricRank"].iloc[relevant_tfs_to_find - 1]
    return r50


print("== BORDA RANKING =============================================")
run_borda()
print("==============================================================")