import sys
import pandas as pd

def readMINI(file, clusters=None):

    """
    convert MINI-EX outupt format to TF-TG edge table
    """

    tmp_tf_list = []
    tmp_tg_list = []
    cluster_list = []

    with open(file, 'r') as inf:
        for i, line in enumerate(inf):
            
            line = line.strip()
            tf = line.split('\t')[0]
            cluster = line.split('\t')[1]
            tgs = set(line.split('\t')[2].split(','))
            
            if clusters != None and cluster not in clusters:
                continue
                
            for tg in tgs:        
                tmp_tf_list.append(tf)
                tmp_tg_list.append(tg)
                cluster_list.append(cluster)

    miniex_df = pd.DataFrame()
    miniex_df['TF'] = tmp_tf_list
    miniex_df['TG'] = tmp_tg_list
    miniex_df['cluster'] = cluster_list
    return miniex_df

def rankMINI(miniex_df, regulons, rank_col='borda_rank'):

    """
    Add regulon borda ranks as edge weights to the TF-TG table
    """

    # prep
    left = miniex_df.set_index(['TF','cluster'])
    right = regulons.set_index(['TF','cluster'])[rank_col].to_frame()

    # join
    return left.join(right, how='left').reset_index()

def scoreMINI(miniex_df, grnboost_df):

    """
    Add GRNBoost2 edge weights to the TF-TG table
    """

    # prep
    left = miniex_df.set_index(['TF','TG'])
    right = grnboost_df.set_index(['TF','TG'])

    # join
    return left.join(right, how='left').reset_index()


if __name__ == '__main__':

    # arguments
    miniex_network_file = sys.argv[1]
    miniex_regulon_file = sys.argv[2]
    grnboost_file = sys.argv[3]
    output_file = sys.argv[4]

    # read GRNBoost2 network
    print("[PROGRESS] reading GRNBoost2 network")
    grnboost_df = pd.read_csv(grnboost_file, sep='\t', header=None)
    grnboost_df.columns = ['TF','TG','weight']

    # read MINI-EX regulon output file
    print("[PROGRESS] reading MINI-EX network")
    regulons = pd.read_excel(miniex_regulon_file)
    regulons['cluster'] = regulons['cluster'].apply(lambda x: '_'.join(x.split('_')[1:]))

    # read MINI-EX network output file
    miniex_df = readMINI(miniex_network_file)

    # add score columns
    print("[PROGRESS] Joining tables")
    miniex_df = rankMINI(miniex_df, regulons, rank_col='borda_rank')
    miniex_df = rankMINI(miniex_df, regulons, rank_col='borda_clusterRank')
    miniex_df = scoreMINI(miniex_df, grnboost_df)

    # sort edges
    print("[PROGRESS] Sorting edges")
    miniex_df['weight_inverse'] = -miniex_df['weight']
    miniex_df = miniex_df.sort_values(['TF','cluster','weight_inverse'])\
                         .drop('weight_inverse', axis=1)

    # write
    print("[PROGRESS] Writing output")
    miniex_df.to_csv(output_file, sep='\t', header=True, index=False)