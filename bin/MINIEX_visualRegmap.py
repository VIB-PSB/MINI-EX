import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.collections as mcoll

from natsort import natsort_keygen

def getCmap(c, n=1000):

    '''
    helper to get fading cmap per tissue type
    '''

    vals = np.ones((n, 4))
    vals[:, 0] = np.linspace(c[0], 1, n)
    vals[:, 1] = np.linspace(c[1], 1, n)
    vals[:, 2] = np.linspace(c[2], 1, n)
    return mcolors.ListedColormap(vals)

def circlemap(data1, data2, group_colors, cluster_grouping, rank_threshold, ax, full_expression=False, alias_info=None):
    
    '''
    Draw heatmap from TF (rows) x Cluster (columns) dataframes
    
    Input:  - data1: cluster rank
            - data2: expression value
            - ax: matplotlib axis to plot on
            - palette
    '''
    
    #circles = []
    for j, cluster in enumerate(data1.columns):
        circles = []
        linewidths = []
        for i, TF in enumerate(data1.index):
            r = data2.loc[TF,cluster]/data2.values.max()/2
            circles.append(plt.Circle((j,i), radius=r))
            linewidths.append(0 if data1.loc[TF,cluster] <= rank_threshold else 0.35)

        if full_expression:
            col = mcoll.PatchCollection(circles, color='white', edgecolor='grey', linewidths=linewidths)
            ax.add_collection(col)
        
        col = mcoll.PatchCollection(circles, match_original=True, cmap=getCmap(group_colors[cluster_grouping[cluster]], n=rank_threshold), linewidths=0.25)
        col.set_array(data1.fillna(rank_threshold+1).loc[:,cluster].values.flatten())
        ax.add_collection(col)

    ax.set_xticks(np.arange(data1.shape[1]))
    ax.set_yticks(np.arange(data1.shape[0]))
    ax.set_xticklabels(data1.columns, rotation=90)
    aliases = [alias_info[g] for g in data1.index.tolist()]
    ax.set_yticklabels(aliases)

    ax.set_xticks((np.arange(data1.shape[1])+1)-0.5, minor=True)
    ax.set_yticks((np.arange(data1.shape[0])+1)-0.5, minor=True)
    ax.grid(which='minor')
    ax.tick_params(axis='both', which='minor', size=0)
    plt.ylim(-0.5, data1.shape[0]-0.5)
    plt.xlim(-0.5, data1.shape[1]-0.5)
    ax.invert_yaxis()
    
def clusterBar(clusters, cluster_types, type_color, ax):
    
    '''
    Plot top row of heatmap: cluster annotation
    
    Input:  - clusters: list of clusters ordered along the trajectory
            - cluster_types: higher level grouping of clusters: dict {<cluster>: <cluster_group>}
            - type_color: dict {<cluster_group>: <color>}
            - ax: matplotlib axis to plot on
    '''
    
    for i, cluster in enumerate(clusters):
        ax.add_patch(plt.Rectangle((i/len(clusters),0), 1/len(clusters), 1, linewidth=0, edgecolor='black', facecolor=type_color[cluster_types[cluster]]))
    ax.set_yticks([])
    ax.set_xticks([])

def getHandles(type_color):
    
    '''
    Helper to build legend
    
    Input:  - type_color: dict {<cluster_group>: <color>}
    Output: - list of handles to pass to plt.legend()
    '''
    
    handles = []
    for cl_type, color in type_color.items():
        handles.append(mpatches.Patch(color=color, label=cl_type))
    return handles

def regMap(regulons, cluster_matrix, clusters, cluster_grouping, group_colors, 
           rank_threshold, outdir, dataset_id, full_expression=False):

    '''
    Plot regulator heatmap

    Input:  - rankedRegulons: pd.DataFrame (MINI-EX output)
            - cluster_matrix: pd.DataFrame (TFs x Clusters)
            - clusters: pd.DataFrame (MINI-EX input: identities, sorted)
            - cluster_grouping: dict {<cluster>: <cluster_group>}
            - group_colors: dict {<cluster_group>: <color>}
            - rank_threshold: int (top x regulators to plot)
            - outdir: str (directory for plots)
            - dataset_id: str (name of dataset, for output name)
            - full_expression: bool (also plot expression for non-top predicted TFs in cluster. default: False)
    Output: - <outdir>/<dataset_id>_regmap_<rank_threshold>.svg
    '''

    # extract alias info
    alias_info = regulons.drop_duplicates('TF').set_index('TF')['alias']
    
    # fetch data
    query = "borda_clusterRank < @rank_threshold"
    data = regulons.query(query)\
                   .pivot_table(index='TF', columns='cluster', values='borda_clusterRank', aggfunc=np.min)\
                   .loc[regulons.query(query)\
                                .sort_values(['cluster','borda_clusterRank'])\
                                .drop_duplicates('TF')\
                                .TF]

    # Init plot
    nrows = len(data)+1
    ncols = len(clusters)
    grid = plt.GridSpec(7, 3, wspace=0, hspace=0, height_ratios=[1,7,3,7,3,7,nrows-28], width_ratios=[ncols,2,1])
    fig = plt.figure(figsize=((ncols+3)/5,nrows/5))

    # plot top cbar
    ax = plt.subplot(grid[0, 0])
    clusterBar(clusters, cluster_grouping, group_colors, ax)
    plt.title('{} | top {} regulators per cluster'.format(dataset_id, rank_threshold))

    # plot heatmap
    ax = plt.subplot(grid[1:, 0])
    circlemap(data, cluster_matrix, group_colors, cluster_grouping, rank_threshold, ax, full_expression=full_expression, alias_info=alias_info)
    ax.set_xticklabels(['{} ({})'.format(cl._text, len(regulons.query("cluster == '{}'".format(cl._text)))) for cl in ax.get_xticklabels()])

    # legend
    ax = plt.subplot(grid[1, 2])
    plt.axis('off')
    ax.scatter([0]*5,[0,1,2,3,4], s=[10,20,40,80,0], c='black')
    plt.text(0.3, 1.5, s='Max expression\nin cluster', rotation=90, ha='center', rotation_mode='anchor')

    ax = plt.subplot(grid[3, 2])
    plt.sca(ax)
    cax = plt.colorbar(cm.ScalarMappable(norm=mcolors.Normalize(vmin=1, vmax=rank_threshold), cmap='Greys_r'), cax=ax)
    cax.set_ticks([1]+[i for i in range(int(rank_threshold/5),rank_threshold+1,int(rank_threshold/5))])
    ax.set_ylabel("Borda rank\n(within cluster)")
    ax.invert_yaxis()

    ax = plt.subplot(grid[5, 2])
    handles = getHandles(group_colors)
    plt.sca(ax)
    plt.legend(handles=handles, bbox_to_anchor=(0,1), loc='upper left', frameon=False, borderpad=0, borderaxespad=0)

    plt.axis('off')
    plt.savefig(f"{outdir}/{dataset_id}_regmap_{rank_threshold}.svg", dpi=200, bbox_inches='tight', facecolor='white')
    # Also provide png for the smaller figures
    if rank_threshold < 40:
        plt.savefig(f"{outdir}/{dataset_id}_regmap_{rank_threshold}.png", dpi=600, bbox_inches='tight', facecolor='white')
    #plt.show()


if __name__ == '__main__':

    # arguments
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--cells2clusters", dest="cells2clusters", type=str, default=None, help="cells2clusters file (MINI-EX INPUT)")
    parser.add_argument("-i", "--identities", dest="identities", type=str, default=None, help="identities file (MINI-EX INPUT)")
    parser.add_argument("-m", "--matrix", dest="matrix", type=str, default=None, help="expression matrix file (MINI-EX INPUT)")
    parser.add_argument("-r", "--regulons", dest="regulons", type=str, default=None, help="regulons file (MINI-EX OUTPUT)")
    parser.add_argument("-t", "--threshold", dest="threshold", type=str, default='25', help="borda_clusterRank threshold (integer or comma-separated integers)")
    parser.add_argument("-p", "--palette", dest="palette", type=str, default='Dark2', help="matplotlib color palette")
    parser.add_argument("-f", "--full_expression", dest="full_expression", type=bool, default=True, help="plot max expression for every cluster")
    parser.add_argument("-o", "--outdir", dest="outdir", type=str, default='.', help="output directory")
    parser.add_argument("-d", "--dataset_id", dest="dataset_id", type=str, default='miniex', help="dataset identifier")
    args = parser.parse_args()

    # check input
    if not args.cells2clusters:
        raise AssertionError("cells2clusters argument must be specified. Run -h flag for help")
    if not args.identities:
        raise AssertionError("identities argument must be specified. Run -h flag for help")
    if not args.matrix:
        raise AssertionError("matrix argument must be specified. Run -h flag for help")
    if not args.regulons:
        raise AssertionError("regulons argument must be specified. Run -h flag for help")

    print("[PROGRESS-REGMAP] Reading inputs")
    sys.stdout.flush()

    # cells2clusters 
    cells2clusters = pd.read_csv(args.cells2clusters, sep='\t', header=None, dtype={0:"str", 1:"str"})\
                       .rename(columns={0:'cell', 1:'cluster'})
    #cells2clusters.columns = ['cell', 'cluster','name']
    cells2clusters.cluster = cells2clusters.cluster.apply(lambda x: 'Cluster_{}'.format(x))
    cells2clusters = cells2clusters.set_index('cell')

    # identities
    identities = pd.read_csv(args.identities, sep='\t', header=None, dtype={0:"str", 1:"str", 2:"int"})
    if identities.shape[1] == 2:
        identities.columns = ['cluster','celltype']
        identities['celltype'] = identities['celltype'].astype(str)
        identities['celltype_lower'] = identities.celltype.str.lower()
        identities = identities.sort_values(['celltype_lower', 'cluster'], key=natsort_keygen())\
                               .drop('celltype_lower', axis=1)
        identities['idx'] = list(range(1,len(identities)+1))
    else:
        identities.columns = ['cluster','celltype','idx']
        identities['celltype'] = identities['celltype'].astype(str)
    identities.cluster = identities.cluster.apply(lambda x: 'Cluster_{}'.format(x))

    # regulons
    regulons = pd.read_excel(args.regulons)
    regulons['cluster'] = regulons['cluster'].apply(lambda x: '_'.join(x.split('_')[-2:]))

    # expression matrix
    cellmatrix = pd.read_csv(args.matrix, sep='\t', index_col=0)
    cellmatrix.index.name = 'geneid'

    # remove non-predicted clusters
    pred_clusters = set(regulons.cluster)
    identities = identities.query("cluster in @pred_clusters")
    cells2clusters = cells2clusters.query("cluster in @pred_clusters")

    # sort clusters
    regulons.cluster = regulons.cluster.astype("category")\
                               .cat.set_categories(identities.sort_values('idx', ascending=True).cluster)

    # collapse expression matrix
    if cellmatrix.shape[1] == len(set(cells2clusters.cluster)):
        cluster_matrix = cellmatrix
    else:
        print("[PROGRESS-REGMAP] Collapsing expression matrix")
        sys.stdout.flush()
        # mean of top three cells per cluster
        cellmatrix = np.log(cellmatrix+1)
        cluster_matrix = cellmatrix.loc[list(set(regulons.TF))]\
                                   .groupby(cells2clusters['cluster'], axis=1)\
                                   .aggregate(lambda v: v.apply(lambda x: x.nlargest(3).mean(), axis=1))
        del cellmatrix

    # define cluster groups
    cluster_grouping = identities.set_index('cluster')['celltype'].to_dict()

    # assign colors
    palette = sns.color_palette(args.palette,len(set(identities.celltype)))
    group_colors = dict(zip(identities.sort_values('idx').drop_duplicates('celltype')['celltype'].tolist(), palette))

    # plot
    threshold_list = sorted(set([int(t) for t in args.threshold.split(',')]+[regulons.borda_clusterRank.max()]))
    for threshold in threshold_list:

        if threshold > regulons.borda_clusterRank.max():
            break

        print("[PROGRESS-REGMAP] Plotting heatmap with top {} regulators".format(threshold))
        sys.stdout.flush()

        regMap(regulons, cluster_matrix, identities.sort_values('idx', ascending=True).cluster, cluster_grouping, 
               group_colors, int(threshold), args.outdir, args.dataset_id, full_expression=args.full_expression)

