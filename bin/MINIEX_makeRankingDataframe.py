"""
Combines selected regulons with all their previously computed metrics into a single dataframe.
"""

import numpy as np
import pandas as pd
import sys

ALIASES_FILE            = sys.argv[1]
CLUSTER_ID_FILE         = sys.argv[2]  
REGULONS_FILE           = sys.argv[3]
CLUSTER_ENRICHMENT_FILE = sys.argv[4]
CENTRALITIES_FILE       = sys.argv[5]
ALL_MARKERS_FILE        = sys.argv[6]
GRNBOOST_FILE           = sys.argv[7]
GO_ENRICHMENT_FILE      = sys.argv[8]
GO_ANNOTATIONS_FILE     = sys.argv[9]
TERMS_OF_INTEREST_FILE  = sys.argv[10]
OUTPUT_FILE             = sys.argv[11]


def path_is_dummy(path): # checks whether the user specified "null" as path to that file
    return path.startswith('.dummy_path')  # file names of unspecified files are replaced with "dummy_path" by the Nextflow pipeline


def retrieve_regulons(ranking_df: pd.DataFrame) -> pd.DataFrame:
    # load regulons filtered according to parameters specified by the user
    regulons_df = pd.read_csv(REGULONS_FILE, sep='\t', names=['TF', 'Cluster', 'TGs'])  # here, clusters are identified as "Cluster_X"
    regulons_df['Regulon'] = regulons_df['TF'] + '_' + regulons_df['Cluster']  # combine TF and cluster names to create a unique regulon name
    regulons_df['#TGs'] = regulons_df['TGs'].str.split(',').str.len()  # count TGs of each regulon

    # collect the relevant columns: 'TF', 'Cluster', '#TGs', 'totRegInCluster'
    ranking_df = regulons_df[['Regulon', 'TF', 'Cluster', 'TGs', '#TGs']].copy()  # collect regulon identifiers + #TGs per regulon
    regulons_per_cluster_dict = regulons_df.groupby('Cluster')['Regulon'].count().to_dict()
    ranking_df['totRegInCluster'] = ranking_df['Cluster'].apply(lambda x: regulons_per_cluster_dict[x])

    return regulons_df


def retrieve_gene_aliases(ranking_df: pd.DataFrame):
    # load gene aliases for the retrieved TFs
    alias_column_names = ['TF', 'alias']
    if not path_is_dummy(ALIASES_FILE):
        tf_alias_df = pd.read_csv(ALIASES_FILE, sep='\t')
        tf_alias_df.columns = alias_column_names
        tf_alias_df = tf_alias_df.fillna("")
    else:
        tf_alias_df = pd.DataFrame(columns=alias_column_names)  # create an empty dataframe
    
    # add entries for the missing TFs, using TF identifiers for the three columns
    missing_tfs = set(ranking_df['TF'].unique()) - set(tf_alias_df['TF'])
    tf_alias_df = pd.concat([tf_alias_df, pd.DataFrame({ 'TF': list(missing_tfs), 'alias': list(missing_tfs)})], ignore_index=True)
    tf_alias_dict = tf_alias_df.groupby('TF')['alias'].agg('/ '.join).to_dict()
    
    # collect the relevant columns: 'alias' and 'clusterId'
    ranking_df['alias'] = ranking_df['TF'].apply(lambda x: tf_alias_dict[x])
    ranking_df['clusterId'] = ranking_df['Cluster'].apply(lambda x: x.split('_')[1])


def retrieve_cluster_annotations(ranking_df: pd.DataFrame):
    # load cluster annotation data
    cluster_id_annotation_df = pd.read_csv(CLUSTER_ID_FILE, sep='\t', names=['clusterId', 'celltype'], dtype={'clusterId': 'str', 'celltype': 'str'})

    # collect the relevant columns: 'cellType' and 'cluster'
    ranking_df = ranking_df.merge(cluster_id_annotation_df, on='clusterId', how='left')
    ranking_df['cluster'] = ranking_df['celltype'] + "_" + ranking_df['Cluster']


def retrieve_cluster_specificity(ranking_df: pd.DataFrame):
    # load cluster specificity for the list of regulons retrieved previously
    cluster_enrichment_df = pd.read_csv(CLUSTER_ENRICHMENT_FILE, sep='\t', comment='#')
    cluster_enrichment_df = cluster_enrichment_df.rename(columns={'set_id': 'TF', 'ftr_id': 'Cluster', 'q-val': 'qval_cluster'})
    cluster_enrichment_df = cluster_enrichment_df[['TF', 'Cluster', 'qval_cluster']]
    
    # collect the relevant columns: 'qval_cluster'
    # only regulons having significant cluster enrichment are kept (inner merge)
    ranking_df = ranking_df.merge(cluster_enrichment_df, on=['TF', 'Cluster'], how='inner')


def retrieve_network_centrality(ranking_df: pd.DataFrame, regulons_df: pd.DataFrame):
    # the original dataframe contains TFs as rows, and, for each cluster, three columns: 'degout_Cluster_X', 'clos_Cluster_X' and 'bet_Cluster_X'
    # plus two columns for the original centrality: 'orig_closeness' and 'orig_betweenness'
    centrality_orig_df = pd.read_csv(CENTRALITIES_FILE, sep='\t', index_col=0)

    # retrieve identifiers (TF + cluster) of the selected regulons
    centrality_df = regulons_df[['TF', 'Cluster']].copy()

    # retrieve original closeness and original betweenness for each TF
    centrality_df = centrality_df.merge(centrality_orig_df, left_on='TF', right_index=True)[['TF', 'Cluster', 'orig_closeness', 'orig_betweenness']]

    # retrieve per cluster information
    centrality_df['out-degree'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"out-degree_{row['Cluster']}"], axis=1)
    centrality_df['cluster_closeness'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"closeness_{row['Cluster']}"], axis=1)
    centrality_df['cluster_betweenness'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"betweenness_{row['Cluster']}"], axis=1)

    # compute final closeness and betweenness:
    # - when defined, use per cluster information
    # - when not defined: use information from the original network
    # - add '+1' to the per cluster information, so that it is always bigger than the original one
    centrality_df["closeness"] = centrality_df["cluster_closeness"] + 1
    centrality_df["betweenness"] = centrality_df["cluster_betweenness"] + 1
    # replace values '1' (meaning: no information available at per cluster level) by information from the original network
    centrality_df["closeness"] = centrality_df["closeness"].replace(1, centrality_df["orig_closeness"])
    centrality_df["betweenness"] = centrality_df["betweenness"].replace(1, centrality_df["orig_betweenness"])

    # collect the final columns
    centrality_df = centrality_df[['TF', 'Cluster', 'out-degree', 'closeness', 'betweenness']]

    # collect the relevant columns: 'out-degree', 'closeness' and 'betweenness'
    ranking_df = ranking_df.merge(centrality_df, on=['TF', 'Cluster'], how='left')


def retrieve_de_information(ranking_df: pd.DataFrame):
    # load cluster DEGs: needed to check whether TFs are DE or not
    markers_df = pd.read_csv(ALL_MARKERS_FILE, sep='\t', dtype={'cluster': str})
    markers_df = markers_df.rename(columns={'gene': 'TF', 'cluster': 'clusterId'})
    markers_df = markers_df[markers_df['p_val_adj'] <= 0.05]

    # collect the relevant columns: 'isTF_DE' and 'TF_qval'
    ranking_df['isTF_DE'] = ranking_df.merge(markers_df, on=['TF', 'clusterId'], how='left', indicator=True)['_merge'].map({'both': 1, 'left_only': 0}).astype(int)
    ranking_df['TF_qval'] = ranking_df.merge(markers_df, on=['TF', 'clusterId'], how='left', indicator=True).apply(lambda row: row['p_val_adj'] if row['_merge'] == 'both' else np.nan, axis=1)


def retrieve_coexpression(ranking_df: pd.DataFrame):
    # load TF-TG co-expression
    grnboost_tf = pd.read_csv(GRNBOOST_FILE, sep='\t', names=['TF', 'TG', 'coexpression'], dtype={'coexpression': float})

    # dictionary associating to each TF a dictionary of its TGs with their coexpression values
    # ex: { 'TF1' : { 'TG1': 0.3, 'TG2': 0.7} ; 'TF2': { 'TG1': 0.2, 'TG4'; 0.1, 'TG5': 0.3 }}
    tf_to_tg_coexpr = (
        grnboost_tf.groupby('TF')
        .apply(lambda x: dict(zip(x['TG'], x['coexpression'])))
        .to_dict()
    )

    # for each row of the regulons dataframe: split its TGs, extracts the corresponding
    # coexpression values and returns the median value
    def compute_median_coexpression(row) -> float:
        tgs = row['TGs'].split(',')
        coexpr_values = [tf_to_tg_coexpr[row['TF']][tg] for tg in tgs]
        return np.median(coexpr_values)
    
    # collect the relevant columns: 'med_coexpr'
    ranking_df['med_coexpr'] = ranking_df.apply(compute_median_coexpression, axis=1)


def retrieve_go_information(ranking_df: pd.DataFrame) -> pd.DataFrame:
    # add default values for the GO-terlated columns
    ranking_df['hasTFrelevantGOterm'] = "unknown_TF"
    ranking_df['GOterm'] = "-"
    ranking_df['GOdescription'] = "-"

    if not path_is_dummy(GO_ENRICHMENT_FILE):  # GO enrichment was performed --> add GO-related information
        # load the list of terms of interest
        if path_is_dummy(TERMS_OF_INTEREST_FILE):
            terms_of_interest = ["DUMMY_TERM_OF_INTEREST"]  # nothing will be matched --> no terms will be considered as relevant
        else:
            terms_of_interest = pd.read_csv(TERMS_OF_INTEREST_FILE, sep='\t', header=None)[0].unique()

        # load GO annotation to determine the list of relevant genes
        go_annotations_df = pd.read_csv(GO_ANNOTATIONS_FILE, sep='\t', header=None, names=['go_term', 'gene_id', 'evidence_code', 'go_description']).sort_values(by=['gene_id', 'go_term'])

        # a term is relevant if either its term or description contain at least one term of interest
        pattern = '|'.join(terms_of_interest)
        go_annotations_df['term_is_relevant'] = go_annotations_df[['go_term', 'go_description']].apply(lambda x: x.str.contains(pattern, case=False, na=False)).any(axis=1)
        all_annotated_genes = set(go_annotations_df['gene_id'])
        annotated_genes_of_interest = set(go_annotations_df[go_annotations_df['term_is_relevant']]['gene_id'])
        ranking_df['hasTFrelevantGOterm'] = ranking_df['TF'].map(
            lambda tf: "relevant_known_TF" if tf in annotated_genes_of_interest
            else "known_TF" if tf in all_annotated_genes
            else "unknown_TF"
        )

        # retrieve GO terms and descriptions
        # for genes associated with at least one relevant term: only relevant terms will be displayed
        # for other genes associated with at least one term: all the terms will be displayed
        def create_go_mapping(df, relevant_only=False) -> dict:
            if relevant_only:
                df = df[df['term_is_relevant']]
            return df.groupby('gene_id').agg(lambda x: ','.join(x)).to_dict()
        
        gene_to_relevant_go_terms_dict = create_go_mapping(go_annotations_df, relevant_only=True)['go_term']
        gene_to_all_go_terms_dict = create_go_mapping(go_annotations_df)['go_term']
        gene_to_relevant_go_descriptions_dict = create_go_mapping(go_annotations_df, relevant_only=True)['go_description']
        gene_to_all_go_descriptions_dict = create_go_mapping(go_annotations_df)['go_description']

        # collect the relevantdata: 'GOterm' and 'GOdescription'
        # add retrieved GO terms and descriptions: first try to associate relevant terms; if not found: all known terms; if not found: '-'
        ranking_df['GOterm'] = ranking_df['TF'].map(gene_to_relevant_go_terms_dict).fillna(ranking_df['TF'].map(gene_to_all_go_terms_dict).fillna('-'))
        ranking_df['GOdescription'] = ranking_df['TF'].map(gene_to_relevant_go_descriptions_dict).fillna(ranking_df['TF'].map(gene_to_all_go_descriptions_dict).fillna('-'))

        # load GO enrichment of target genes
        go_enrichment_df = pd.read_csv(GO_ENRICHMENT_FILE, sep='\t', comment='#')
        go_enrichment_df = go_enrichment_df.rename(columns={'set_id': 'Regulon', 'ftr_id': 'GO_enrich_term', 'q-val': 'GO_enrich_qval'})
        go_enrichment_df = go_enrichment_df[go_enrichment_df['GO_enrich_term'].isin(set(go_annotations_df[go_annotations_df['term_is_relevant']]['go_term']))]  # only select relevant GO terms
        go_enrichment_df['#TGs_withGO'] = go_enrichment_df['hits'].str.split(',').str.len()  # how many TGs are associated with the GO term
        go_enrichment_df = go_enrichment_df.merge(go_annotations_df[['go_term', 'go_description']], left_on='GO_enrich_term', right_on='go_term')  # add GO description
        go_enrichment_df = go_enrichment_df.rename(columns={'go_description': 'GO_enrich_desc'})
        go_enrichment_df = go_enrichment_df[['Regulon', 'GO_enrich_qval', 'GO_enrich_term', 'GO_enrich_desc', '#TGs_withGO']]  # select relevant columns
        go_enrichment_df = go_enrichment_df.loc[go_enrichment_df.groupby('Regulon')['GO_enrich_qval'].idxmin()]  # select elements with the lowest enrichment q-value
        
        # collect the relevant data: 'GO_enrich_qval', 'GO_enrich_term', 'GO_enrich_desc' and '#TGs_withGO'
        ranking_df = ranking_df.merge(go_enrichment_df, on=['Regulon'], how='left')
        ranking_df[['GO_enrich_term', 'GO_enrich_desc']] = ranking_df[['GO_enrich_term', 'GO_enrich_desc']].fillna('-')


    # collect the relevant columns columns and reorder them
    if path_is_dummy(TERMS_OF_INTEREST_FILE):  # ignore GO-related columns
        ranking_df = ranking_df[['TF', 'alias', 'hasTFrelevantGOterm', 'GOterm', 'GOdescription', 'cluster', 'celltype', 'isTF_DE',
                                'totRegInCluster', '#TGs', 'qval_cluster', 'out-degree', 'closeness', 'betweenness', 'med_coexpr', 'TF_qval']]
    else:  # add GO-related columns
        ranking_df = ranking_df[['TF', 'alias', 'hasTFrelevantGOterm', 'GOterm', 'GOdescription', 'cluster', 'celltype', 'isTF_DE', 'totRegInCluster', '#TGs',
                                'qval_cluster', 'out-degree', 'closeness', 'betweenness', 'med_coexpr', 'TF_qval', 'GO_enrich_qval', 'GO_enrich_term', 'GO_enrich_desc', '#TGs_withGO']]
        
    return go_annotations_df


# ======== dataframe that will collect all the metrics computed for the selected regulons ========
ranking_df = pd.DataFrame()

# ======== collect the information from different sources ========
regulons_df = retrieve_regulons(ranking_df)
retrieve_gene_aliases(ranking_df)
retrieve_cluster_annotations(ranking_df)
retrieve_cluster_specificity(ranking_df)
retrieve_network_centrality(ranking_df, regulons_df)
retrieve_de_information(ranking_df)
retrieve_coexpression(ranking_df)
go_annotations_df = retrieve_go_information(ranking_df)

# ======== save the resulting dataframe ========
ranking_df.to_csv(OUTPUT_FILE, sep='\t', index=None)

# ======== write log information (only if the user provided terms of interest) ========
if not path_is_dummy(TERMS_OF_INTEREST_FILE):
    print("== SELECTED GO TERMS =========================================")
    print("Following GO terms were retrieved for selected terms of interest:")
    with pd.option_context('display.max_rows', None):  # do not truncate rows for large dataframes
        print(go_annotations_df[go_annotations_df['term_is_relevant'] == True][['go_term', 'go_description']].drop_duplicates().sort_values(by='go_term').to_string(index=False))
    print("")