import numpy,pandas,collections,statistics,natsort,sys


TF_MATRIX=sys.argv[1]
TF_NAMES=sys.argv[2]
CELL_NAMES=sys.argv[3]
GRNBOOST_GRN=sys.argv[4]
MOTIF_ENRICHMENT_GRN=sys.argv[5]
EXPRESSION_FILTERED_GRN=sys.argv[6]
TF_FILE=sys.argv[7]
CELLS=sys.argv[8]
CLUSTER_IDS=sys.argv[9]
OUT=sys.argv[10]

####################
### READ IN DATA ###
####################

# Read in TFs as set
tfs_all = set(pandas.read_csv(TF_FILE, header=None, sep='\t')[0])

# Read in count matrix
with open(TF_NAMES, 'r') as f:
    gene_names = [line.strip() for line in f]
with open(CELL_NAMES, 'r') as f:
    cell_names = [line.strip() for line in f]
count_matrix = pandas.DataFrame(numpy.load(TF_MATRIX).T, index=gene_names, columns=cell_names)

# Read in GRNBoost2 output as TF-to-nb_of_TGs dict
tf2nb_of_tgs_grnboost = collections.Counter(pandas.read_csv(GRNBOOST_GRN, header=None, sep='\t')[0])

# Read in motif enrichment output as TF-to-nb_of_TGs dict
tf2nb_of_tgs_motif_enrichment = collections.Counter(pandas.read_csv(MOTIF_ENRICHMENT_GRN, header=None, sep='\t')[0])

# Read in expression filtered output (final regulons) as a dict with cluster as key and a TF-to-nb_of_TGs dict as value
cluster2tf2nb_of_tgs_expression_filtered = dict()
with open(EXPRESSION_FILTERED_GRN) as file:
    for line in file:
        tf, cluster, tgs = line.rstrip().split('\t')
        if cluster not in cluster2tf2nb_of_tgs_expression_filtered:
            cluster2tf2nb_of_tgs_expression_filtered[cluster] = dict()
        # add the TF as key to the dict, with the number of TGs (as nb of comma's + 1) as value
        cluster2tf2nb_of_tgs_expression_filtered[cluster][tf] = tgs.count(',') + 1

# Read in cell2cluster info
cell_cluster_df = pandas.read_csv(CELLS, sep='\t', header=None, names=["cell_id", "cluster"], usecols=["cell_id", "cluster"], index_col="cell_id", dtype={"cluster":"str"})

# Read in cluster (number) and tissue identity and merge into one string (e.g. "xylem-28")
cluster_ids_df = pandas.read_csv(CLUSTER_IDS, sep='\t', header=None, names=["cluster_id", "tissue"], dtype={"cluster_id":"str", "tissue":"str"})
cluster_ids_df["merged_cluster_name"] = cluster_ids_df['tissue'].astype(str) + '-' + cluster_ids_df['cluster_id'].astype(str)
cluster_ids_df.set_index("cluster_id", inplace=True)

################################
### CREATE FILTERING INFO DF ###
################################

# Get all expressed TFs
tfs_expressed = count_matrix.index

# Get all TFs of the final expression filtered regulons
tfs_expression_filtered = set().union(*([set(tf2nb.keys()) for tf2nb in cluster2tf2nb_of_tgs_expression_filtered.values()]))

# Make df with info on which TFs are kept at each filtering step
# Start with a 0-filled df
filtering_steps = ['isTF_expressed','isTF_expressionRegulons','isTF_TFBSRegulons','isTF_finalRegulons']
filtering_info_df = pandas.DataFrame(0, index=sorted(tfs_all), columns=filtering_steps)

# Fill the df with 1's if the TF is not filtered out at this step
for tf in filtering_info_df.index:
    if tf in tfs_expressed:
        filtering_info_df.loc[tf,'isTF_expressed'] = 1
    if tf in tf2nb_of_tgs_grnboost.keys():
        filtering_info_df.loc[tf,'isTF_expressionRegulons'] = 1
    if tf in tf2nb_of_tgs_motif_enrichment.keys():
        filtering_info_df.loc[tf,'isTF_TFBSRegulons'] = 1
    if tf in tfs_expression_filtered:
        filtering_info_df.loc[tf,'isTF_finalRegulons'] = 1

#################################
### CREATE EXPRESSION INFO DF ###
#################################

# Subset TF matrix to only keep annotated cells
count_matrix = count_matrix.transpose()
count_matrix = count_matrix[count_matrix.index.isin(cell_cluster_df.index)]

# Get a dict with cluster ID as keys and number of cells as values
cluster2nb_of_cells = collections.Counter(cell_cluster_df.cluster)

# Start with a 0-filled df 
expression_info_df = pandas.DataFrame(0.0, index=sorted(tfs_expressed), columns=sorted(cluster_ids_df["merged_cluster_name"]))

for tf in tfs_expressed:
    # Make a dict with cluster ID as keys and number of cells expressing the TF as values
    cells_expressing_tf = count_matrix[count_matrix[tf] > 1].index
    cluster2nb_of_cells_with_tf = collections.Counter(cell_cluster_df.loc[cells_expressing_tf].cluster)

    # For each cluster, calculate the percentage of cells expressing the TF
    for cluster, nb_of_cells_with_tf in cluster2nb_of_cells_with_tf.items():
        percentage = round((nb_of_cells_with_tf / cluster2nb_of_cells[cluster]) * 100, 1)
        # Add the percentage to the expression info df
        expression_info_df.loc[tf, cluster_ids_df.loc[cluster, "merged_cluster_name"]] = percentage

################################
### COMBINE AND WRITE OUTPUT ###
################################

# Add the expression info to the filtering info df
complete_info_df = filtering_info_df.join(expression_info_df)
# Write as output file
complete_info_df.to_csv(OUT, sep='\t')

#################################
### EXTRACT INFO FOR LOG FILE ###
#################################

# Display the number of initial TFs and the number of expressed TFs
data = {'Number of TFs': [len(tfs_all), len(tfs_expressed)]}
print("== INITIAL TF SET ============================================")
print(pandas.DataFrame(data, index=['All TFs', 'Expressed TFs']))
print("")

# Get all regulon numbers and regulon sizes for all clusters
cluster2nb_of_regulons = {cluster: len(tf2nb_of_tgs) for cluster, tf2nb_of_tgs in cluster2tf2nb_of_tgs_expression_filtered.items()}
target_gene_counts = [count for cluster_dict in cluster2tf2nb_of_tgs_expression_filtered.values() for count in cluster_dict.values()]
cluster2median_size = {cluster: statistics.median(tf2nb_of_tgs.values()) for cluster, tf2nb_of_tgs in cluster2tf2nb_of_tgs_expression_filtered.items()}

# Prepare data to display in the log file for each filtering step
data = {'Unique regulators': [len(set(tf2nb_of_tgs_grnboost.keys())), \
                              len(set(tf2nb_of_tgs_motif_enrichment.keys())), \
                              len(tfs_expression_filtered)],
        'Regulons': [len(tf2nb_of_tgs_grnboost), \
                     len(tf2nb_of_tgs_motif_enrichment), \
                     sum(cluster2nb_of_regulons.values())],
        'Median regulon size': [statistics.median(tf2nb_of_tgs_grnboost.values()), \
                                statistics.median(tf2nb_of_tgs_motif_enrichment.values()), \
                                statistics.median(target_gene_counts)]}

# Update Pandas display parameters so that the dataframes are not truncated
pandas.set_option('display.max_columns', None)
pandas.set_option('display.width', 1000)
print("== GRN FILTERING =============================================")
print(pandas.DataFrame(data, index=['Step 1 (GRNBoost2)', 'Step 2 (motif filtering)', 'Step 3 (expression filtering)']))
print("")
print("Per cluster statistics of the final GRN:")
print("")
clusters = natsort.natsorted(cluster2nb_of_regulons.keys())
data = {'Number of regulons': [cluster2nb_of_regulons[cluster] for cluster in clusters], \
        'Regulon size (median number of TGs)': [cluster2median_size[cluster] for cluster in clusters]}
per_cluster_df = pandas.DataFrame(data, index=clusters)
medians = per_cluster_df.median()
per_cluster_df = per_cluster_df.append(medians.rename('MEDIAN VALUE'))
print(per_cluster_df)
print("")
