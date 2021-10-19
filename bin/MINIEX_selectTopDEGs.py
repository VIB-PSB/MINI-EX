import pandas,sys

CLUSTERS=sys.argv[1]
df=pandas.read_csv(CLUSTERS,sep='\t',dtype={'cluster': str})
df=df[df['p_val_adj']< 0.05]

def select_topDEGs(t,o):       
    t=int(t)
    df_top=df.sort_values('p_val_adj',ascending=True).groupby('cluster').head(t) 

    fout = df_top[['cluster', 'gene']].copy()
    fout['cluster'] = 'Cluster_'+fout['cluster']
    fout.to_csv(o,sep='\t',header=None,index=None)

tops=sys.argv[2]
outfile=sys.argv[3]
select_topDEGs(tops, outfile)
