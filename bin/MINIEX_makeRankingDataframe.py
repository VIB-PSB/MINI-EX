"""
Combines selected regulons with all their previously computed metrics into a single dataframe.
"""

import pandas as pd
import sys

ALIASES_FILE            = sys.argv[1]
CLUSTER_ID_FILE         = sys.argv[2]  
REGULONS_FILE           = sys.argv[3]
CLUSTER_ENRICHMENT_FILE = sys.argv[4]
CENTRALITIES_FILE       = sys.argv[5]
ALL_MARKERS_FILE        = sys.argv[6]
GO_ENRICHMENT_FILE      = sys.argv[7]
GO_ANNOTATIONS_FILE     = sys.argv[8]
TERMS_OF_INTEREST_FILE  = sys.argv[9]
OUTPUT_FILE             = sys.argv[10]


def path_is_dummy(path): # checks whether the user specified "null" as path to that file
    return path.startswith('.dummy_path')  # file names of unspecified files are replaced with "dummy_path" by the Nextflow pipeline


# ======== dataframe that will collect all the metrics computed for the selected regulons ========
ranking_df = pd.DataFrame()


# ======== load regulons filtered according to parameters specified by the user ========
regulons_df = pd.read_csv(REGULONS_FILE, sep='\t', names=['TF', 'Cluster', 'TGs'])  # here, clusters are identified as "Cluster_X"
regulons_df['Regulon'] = regulons_df['TF'] + '_' + regulons_df['Cluster']  # combine TF and cluster names to create a unique regulon name
regulons_df['#TGs'] = regulons_df['TGs'].str.split(',').str.len()  # count TGs of each regulon
regulons_dict = regulons_df.set_index('Regulon')['#TGs'].to_dict()
# collect the relevant data: 'TF', 'Cluster', '#TGs', 'totRegInCluster'
ranking_df = regulons_df[['Regulon', 'TF', 'Cluster', '#TGs']].copy()  # collect regulon identifiers + #TGs per regulon
regulons_per_cluster_dict = regulons_df.groupby('Cluster')['Regulon'].count().to_dict()
ranking_df['totRegInCluster'] = ranking_df['Cluster'].apply(lambda x: regulons_per_cluster_dict[x])


# ======== load gene aliases for the retrieved TFs ========
alias_column_names = ['TF', 'alias', 'description']
if not path_is_dummy(ALIASES_FILE):
    tf_alias_df = pd.read_csv(ALIASES_FILE, sep='\t')  # original columns: 'locus_name' - 'symbol' - 'full_name'
    tf_alias_df.columns = alias_column_names
    tf_alias_df = tf_alias_df.fillna("")
else:
    tf_alias_df = pd.DataFrame(columns=alias_column_names)  # create an empty dataframe
# add entries for the missing TFs, using TF identifiers for the three columns
missing_tfs = set(ranking_df['TF'].unique()) - set(tf_alias_df['TF'])
tf_alias_df = pd.concat([tf_alias_df, pd.DataFrame({ 'TF': list(missing_tfs), 'alias': list(missing_tfs), 'description': list(missing_tfs)})], ignore_index=True)
tf_alias_dict = tf_alias_df.groupby('TF')['alias'].agg('/ '.join).to_dict()
# collect the relevant data: 'alias' and 'clusterId'
ranking_df['alias'] = ranking_df['TF'].apply(lambda x: tf_alias_dict[x])
ranking_df['clusterId'] = ranking_df['Cluster'].apply(lambda x: x.split('_')[1])


# ======== load cluster annotations ========
cluster_id_annotation_df = pd.read_csv(CLUSTER_ID_FILE, sep='\t', names=['clusterId', 'celltype'], dtype={'clusterId': 'str'})
# collect the relevant data: 'cellType' and 'cluster'
ranking_df = ranking_df.merge(cluster_id_annotation_df, on='clusterId', how='left')
ranking_df['cluster'] = ranking_df['celltype'] + "_" + ranking_df['Cluster']


# ======== load cluster specificity for the list of regulons retrieved previously ========
cluster_enrichment_df = pd.read_csv(CLUSTER_ENRICHMENT_FILE, sep='\t', comment='#')
cluster_enrichment_df = cluster_enrichment_df.rename(columns={'set_id': 'TF', 'ftr_id': 'Cluster', 'q-val': 'qval_cluster'})
cluster_enrichment_df = cluster_enrichment_df[['TF', 'Cluster', 'qval_cluster']]
# collect the relevant data: 'qval_cluster'
# only regulons having significant cluster enrichment are kept (inner merge)
ranking_df = ranking_df.merge(cluster_enrichment_df, on=['TF', 'Cluster'], how='inner')


# ======== load network centrality measures ========
# the original dataframe contains TFs as rows, and, for each cluster, three columns: 'degout_Cluster_X', 'clos_Cluster_X' and 'bet_Cluster_X'
centrality_orig_df = pd.read_csv(CENTRALITIES_FILE, sep='\t', index_col=0)
centrality_df = regulons_df[['TF', 'Cluster']].copy()  # retrieve identifiers of the selected regulons
centrality_df['out-degree'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"degout_{row['Cluster']}"], axis=1)
centrality_df['closeness'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"clos_{row['Cluster']}"], axis=1)
centrality_df['betweenness'] = centrality_df.apply(lambda row: centrality_orig_df.loc[row['TF'], f"bet_{row['Cluster']}"], axis=1)
# collect the relevant data: 'out-degree', 'closeness' and 'betweenness'
ranking_df = ranking_df.merge(centrality_df, on=['TF', 'Cluster'], how='left')


# ======== load cluster DEGs ========
# this is needed to check whether TFs are DE or not
markers_df = pd.read_csv(ALL_MARKERS_FILE, sep='\t',dtype={'cluster': str})
markers_df = markers_df.rename(columns={'gene': 'TF', 'cluster': 'clusterId'})
markers_df = markers_df[markers_df['p_val_adj'] <= 0.05]
# collect the relevant data: 'isTF_DE'
ranking_df['isTF_DE'] = ranking_df.merge(markers_df, on=['TF', 'clusterId'], how='left', indicator=True)['_merge'].map({'both': 1, 'left_only': 0}).astype(int)


# ======== add GO-related information ========
# add default values for the GO-terlated columns
ranking_df['hasTFrelevantGOterm'] = "unknown_TF"
ranking_df['GOterm'] = "-"
ranking_df['GOdescription'] = "-"


if path_is_dummy(GO_ENRICHMENT_FILE):  # no GO enrichment was performed --> keep the default values
    # select final columns and reorder them
    ranking_df = ranking_df[['TF', 'alias', 'hasTFrelevantGOterm', 'GOterm', 'GOdescription', 'cluster', 'celltype', 'isTF_DE',
                             'totRegInCluster', '#TGs', 'qval_cluster', 'out-degree', 'closeness', 'betweenness']]

else:  # add GO-related columns
    # ======== load the list of terms of interest ========
    if path_is_dummy(TERMS_OF_INTEREST_FILE):
        terms_of_interest = ["DUMMY_TERM_OF_INTEREST"]  # nothing will be matched --> no terms will be considered as relevant
    else:
        terms_of_interest = pd.read_csv(TERMS_OF_INTEREST_FILE, sep='\t', header=None)[0].unique()

    # ======== load GO annotation to determine the list of relevant genes ========
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
    def create_go_mapping(df, relevant_only=False):
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

    # ======== load GO enrichment of target genes ========
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

    # select final columns and reorder them
    ranking_df = ranking_df[['TF', 'alias', 'hasTFrelevantGOterm', 'GOterm', 'GOdescription', 'cluster', 'celltype', 'isTF_DE', 'totRegInCluster', '#TGs',
                             'qval_cluster', 'out-degree', 'closeness', 'betweenness', 'GO_enrich_qval', 'GO_enrich_term', 'GO_enrich_desc', '#TGs_withGO']]


# ======== save the resulting dataframe ========
ranking_df.to_csv(OUTPUT_FILE, sep='\t', index=None)


# write log information (only if the user provided terms of interest)
if not path_is_dummy(TERMS_OF_INTEREST_FILE):
    print("== SELECTED GO TERMS =========================================")
    print("Following GO terms were retrieved for selected terms of interest:")
    with pd.option_context('display.max_rows', None):  # do not truncate rows for large dataframes
        print(go_annotations_df[go_annotations_df['term_is_relevant'] == True][['go_term', 'go_description']].drop_duplicates().sort_values(by='go_term').to_string(index=False))
    print("")