from itertools import combinations
from typing import Dict
from scipy import stats
import sys
import pandas as pd

REGULONS_FILE = sys.argv[1]
OUTPUT_FILE   = sys.argv[2]
PROCEDURE     = sys.argv[3]


def run_borda():
    # the regulons dataframe provides all the information about the regulons
    # collected so far, including the values of the individual metrics
    regulons_df = pd.read_csv(REGULONS_FILE, sep='\t')

    # list of metrics that can be used for borda ranking + their sorting order ('ascending' = lower values have better rank)
    metrics = {'qval_cluster': 'ascending', 'out-degree': 'descending', 'betweenness': 'descending', 'closeness': 'descending', 'TF_qval': 'ascending', 'med_coexpr': 'descending', 'GO_enrich_qval': 'ascending'}
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

    regulons_df.to_excel(f"{OUTPUT_FILE}.xlsx", index=None)
    regulons_df.to_csv(f"{OUTPUT_FILE}.tsv", sep='\t', index=None)


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
    - for each metric, computes its scaling factor (=1 for the metrics with the lowest R50,
      =0 for the metric assigning random ranks)
    - discards metrics with scaling factors <= 0
    - for each retained metric, computes its weighted ranks by multiplying the original
      ranks by their scaling factor
    - for each combination of retained metrics (going from a single metric to all the metrics),
      computes the sum of the weighted ranks of the individual metrics, sorts them from lowest to highest,
      and computes the R50 corresponding to that ranking (=Borda R50)
    - select the best combination of metrics, defined as having the lowest Borda R50
    - using the selected combination of metrics, computes the global Borda rank (resulting
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

    # get R50 values for each metric
    metric_r50_dict = { metric: get_r50_for_metric(regulons_df, metric) for metric in metrics.keys() }

    # identify best (= lowest R50) and worst (half of unique regulons) ranks
    best_rank = min(metric_r50_dict.values())
    worst_rank = round(regulons_df['TF'].nunique() / 2)

    # compute scaling factor for each metric (=1 for the best metric, =0 for the worst metric, and between 0 and 1 for everything else)
    metric_r50_scaling_factor_dict = {key: max(0, (worst_rank - rank) / (worst_rank - best_rank)) for key, rank in metric_r50_dict.items() }

    # compute weighted ranks, by multiplying the ranks by the scaling factor
    for metric, scaling_factor in metric_r50_scaling_factor_dict.items():
        regulons_df[f'{metric}_metricRankWeighted'] = regulons_df[f'{metric}_metricRank'] * scaling_factor
        regulons_df[f'{metric}_metricClusterRankWeighted'] = regulons_df[f'{metric}_metricClusterRank'] * scaling_factor

    # only select metrics with positive scaling factors
    selected_metrics = [metric for metric, scaling_factor in metric_r50_scaling_factor_dict.items() if scaling_factor > 0]

    # generate all the combinations of the metrics
    metric_combinations = [comb for i in range(1, len(selected_metrics) + 1) for comb in combinations(selected_metrics, i)]

    # print R50 in the log
    log_df = pd.DataFrame(columns=['R50'])
    for metric, r50 in metric_r50_dict.items():
        log_df.loc[metric] = f"{r50}"
    log_df.loc["----------------------------"] = "---"
    log_df.loc["best R50"] = f"{best_rank}"
    log_df.loc["worst R50 (= random ranking)"] = f"{worst_rank}"
    print(log_df)

    # compute Borda R50 for each combination of metrics and store it into a dictionary
    metric_combination_r50 = {}
    for metric_combination in metric_combinations:
        weighted_metric_combination = [item + '_metricRankWeighted' for item in metric_combination]
        regulons_df['borda'] = regulons_df[weighted_metric_combination].sum(axis=1)
        regulons_df['borda_metricRank']=regulons_df['borda'].rank(method='max', ascending=True)
        borda_r50 = get_r50_for_metric(regulons_df, 'borda')
        metric_combination_r50[','.join(metric_combination)] = borda_r50    

    # compute global Borda rank
    best_metric_combination_found = False
    while not best_metric_combination_found:
         # identify the best combination of metrics (=having the lowest Borda R50)
        best_metric_combination = min(metric_combination_r50, key=lambda k: metric_combination_r50[k])
        best_metric_combination_list = best_metric_combination.split(',')

        # compute final Borda rank based on the best combination of metrics
        weighted_metric_combination = [item + '_metricRankWeighted' for item in best_metric_combination_list]
        regulons_df['borda'] = regulons_df[weighted_metric_combination].sum(axis=1)
        regulons_df['borda_rank'] = regulons_df['borda'].rank(method='max', ascending=True)

        # ranks are valid if the first 10 or the last 50% of the regulons are not assigned to the same Borda rank
        if borda_ranks_are_valid(regulons_df, best_metric_combination):
            best_metric_combination_found = True
            # add the information about the selected metrics to the log file
            print(f"R50 for the selected metrics  {metric_combination_r50[best_metric_combination]}")
            print("----------------------------  ---")
            print(f"Selected metrics              {best_metric_combination}")
        else:
            # remove that combination of the dictionary and check the next best combination
            del metric_combination_r50[best_metric_combination]
    
    # compute final Borda rank per cluster based on the best combination of metrics
    weighted_cluster_metric_combination = [item + '_metricClusterRankWeighted' for item in best_metric_combination_list]
    regulons_df['bordaCluster'] = regulons_df[weighted_cluster_metric_combination].sum(axis=1)
    regulons_df['borda_clusterRank'] = regulons_df.groupby('cluster')['bordaCluster'].rank(method='max', ascending=True)

    return regulons_df


def borda_ranks_are_valid(regulons_df: pd.DataFrame, metric_combination: str) -> bool:
    """
    Validates the Borda ranking for the provided metric combination.

    Borda ranking can be invalid if:
    - the first 10 elements are assigned to the same rank
    - the last 50% of elements are assigned to the same rank

    Parameters
    ----------
    regulons_df : pd.DataFrame
        The dataframe containing the regulons and their Borda ranks.
    metric_combination : str
        The metric combination being evaluated.

    Returns
    -------
    bool
        True if the Borda ranks are valid, False if there are ties in the first
        10 ranks or the last 50% of ranks.
    """
    sorted_df = regulons_df.sort_values(by="borda_rank")
    
    first_elements = sorted_df['borda_rank'].iloc[:10]
    if first_elements.nunique() == 1:
        print(f"WARNING: the best metric combination '{metric_combination}' was discarded as it assigned equal ranks for the first 10 regulons.")
        return False
    
    last_elements = sorted_df['borda_rank'].iloc[(len(sorted_df) // 2):]
    if last_elements.nunique() == 1:
        print(f"WARNING: the best metric combination '{metric_combination}' was discarded as it assigned equal ranks for the last 50% of regulons.")
        return False
    
    return True


def get_r50_for_metric(regulons_df: pd.DataFrame, metric_name: str):
    """
    Computes the R50 value for the provided metric.
    R50 corresponds to the rank at which half of the TFs
    annotated as relevant are recovered using the provided metric.
    Only unique TFs are considered: if a TF is present in several clusters,
    its position corresponding to the lowest rank is retrieved.
    The selected TFs then received reassigned ranks: from 1 to the number of unique TFs.

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
    # for each TF: get indices corresponding to their lowest ranks and retrieve the list of selected TFs
    min_rank_indices = regulons_df.groupby("TF")[f"{metric_name}_metricRank"].idxmin()
    unique_tfs_df = regulons_df.loc[min_rank_indices]

    # sort the selected TFs by the original metric rank
    sorted_unique_tfs_df = unique_tfs_df.sort_values(by=f"{metric_name}_metricRank")

    # reassign ranks incrementally (1-based indexing)
    sorted_unique_tfs_df['ReassignedRank'] = range(1, len(sorted_unique_tfs_df) + 1)

    # select the relevant TFs
    relevant_tfs_df = sorted_unique_tfs_df[sorted_unique_tfs_df['hasTFrelevantGOterm'] == 'relevant_known_TF']

    # compute the median reassigned rank for the relevant TFs
    relevant_tfs_to_find = round(len(relevant_tfs_df) / 2)
    r50 = relevant_tfs_df['ReassignedRank'].iloc[relevant_tfs_to_find - 1]

    return r50


print("== BORDA RANKING =============================================")
run_borda()
print("")
print("==============================================================")