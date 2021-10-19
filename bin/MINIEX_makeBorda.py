import pandas, sys, os
from itertools import combinations
from scipy import stats


def make_ranks(lisOfTerms,df_start):
    tfclu=[row['TF']+'_'+row['cluster'] for index,row in df_start.iterrows()]
    ranks_df=pandas.DataFrame(index=tfclu)
    for w in lisOfTerms:
        sing_rank={}
        for index,row in df_start.iterrows():
            sing_rank[row['TF']+'_'+row['cluster']]=row[w]
        sing_rank=pandas.DataFrame.from_dict(sing_rank, orient='index')
        if w == 'qval_cluster' or w == 'GO_enrich_qval': 
            sing_rank[w]=sing_rank[0].rank(method ='max',ascending=True,na_option='bottom')
        else:
            sing_rank[w]=sing_rank[0].rank(method ='max',ascending=False,na_option='bottom')
        sing_rank=sing_rank.drop([0], axis=1)
        ranks_df=ranks_df.join(sing_rank, how='left')    
    return ranks_df


def cutoff(metric,ranked_df):
    ranked_df=ranked_df.sort_values(by=metric)
    trueObs=0
    for index,row in ranked_df.iterrows():
        if '_'.join(index.rsplit('_')[:-3]) in myTFs:
            trueObs+=1
            if trueObs == ct:
                return row[metric]




def calculate_borda_ref(df_in,df_ref):
    df_out = pandas.concat([pandas.Series(df_in[k]/cutoff(k,df_ref)) for k in list(df_in.columns)], axis=1)
    df_out['borda'] = df_out.iloc[: ,:].sum(axis=1) 
    df_out['borda_rank']=df_out['borda'].rank(method ='max',ascending=True)
    df_out=df_out.drop(['borda'], axis=1)
    df_out=df_out.sort_values(by='borda_rank')   
    return df_out

def borda_cluster_ref(df_in,combina,ref_df):    
    clu_ranks=[]
    for clu in list(set(df_in['cluster'])):
        tmp_rank=calculate_borda_ref(make_ranks(combina,df_in[df_in['cluster']==clu]),ref_df)
        clu_ranks.append(tmp_rank['borda_rank'])
    return pandas.concat(clu_ranks)    

def calculate_borda_std(df_in):
    df_out=pandas.DataFrame(index=df_in.index)
    df_out['borda'] = stats.gmean(df_in.iloc[: ,:], axis=1)
    df_out['borda_rank']=df_out['borda'].rank(method ='max',ascending=True)
    df_out=df_out.drop(['borda'], axis=1)
    df_out=df_out.sort_values(by='borda_rank')   
    return df_out

def borda_cluster_std(df_in):    
    clu_ranks=[]
    for clu in list(set(df_in['cluster'])):
        tmp_rank=calculate_borda_std(make_ranks(METRICS,df_in[df_in['cluster']==clu]))
        clu_ranks.append(tmp_rank['borda_rank']) 
    return pandas.concat(clu_ranks)     




DF_IN=sys.argv[1]
METRICS=sys.argv[2]
DF_OUT=sys.argv[3]
REF=sys.argv[4]

df=pandas.read_csv(DF_IN,sep='\t')
METRICS=METRICS.rstrip().rsplit(',')

if REF == 'ref':
    
    myTFs=list(set(df[(df['hasTFrelevantGOterm']=='relevant_known_TF')]["TF"]))
    ct=round(list(df['hasTFrelevantGOterm']).count('relevant_known_TF')/2)
   
    
    ###for each metric do the ranking (same vals get subsequential ranks)        
    combis = sum([list(map(list, combinations(METRICS, i))) for i in range(len(METRICS) + 1)], [])
    tmp_dic={}
    for c in combis:
        if len(c)>0:
            df_rank=make_ranks(c,df)
            df_wei=calculate_borda_ref(df_rank,df_rank)
            
            cuto=cutoff('borda_rank',df_wei)
            tmp_dic[','.join(c)]=[c,cuto,df_wei['borda_rank'],df_rank]
    
    
    
    best_combi=min(list(tmp_dic.values()), key=lambda x: x[1])[0]
    print('Calculating BORDA ranks for %s using the best metrics combination: '% os.path.basename(DF_IN).rsplit('_')[0], ', '.join(best_combi))
    ranks=pandas.DataFrame(min(list(tmp_dic.values()), key=lambda x: x[1])[2]).copy()
    reference_df=pandas.DataFrame(min(list(tmp_dic.values()), key=lambda x: x[1])[3]).copy()
    
    ranks['index1']=ranks.index
    ranks['TF']=ranks.index1.str.split('_').str[:-3].str.join("_")
    ranks['cluster']=ranks.index1.str.split('_').str[-3:].str.join("_")
    ranks=ranks.drop(['index1'],axis=1)
    
    mid=pandas.DataFrame(borda_cluster_ref(df,best_combi,reference_df))
    mid.columns=['borda_clusterRank']
    
    ranks['borda_clusterRank'] = mid
    df_final=pandas.merge(df,ranks[["TF", "cluster","borda_rank","borda_clusterRank"]],on=["TF", "cluster"], how='left')
    
    df_final=df_final.sort_values(by="borda_rank")
    
    df_final.to_excel(DF_OUT,index=None)

elif REF == 'std':
    print('Calculating BORDA ranks for %s on default parameters including : '% os.path.basename(DF_IN).rsplit('_')[0], ', '.join(METRICS))
    ranks=calculate_borda_std(make_ranks(METRICS,df))
        
    ranks['index1']=ranks.index
    ranks['TF']=ranks.index1.str.split('_').str[:-3].str.join("_")
    ranks['cluster']=ranks.index1.str.split('_').str[-3:].str.join("_")
    ranks=ranks.drop(['index1'],axis=1)
    
    mid=pandas.DataFrame(borda_cluster_std(df))
    mid.columns=['borda_clusterRank']
    
    ranks['borda_clusterRank'] = mid
    df_final=pandas.merge(df,ranks[["TF", "cluster","borda_rank","borda_clusterRank"]],on=["TF", "cluster"], how='left')
    
    df_final=df_final.sort_values(by="borda_rank")
    
    df_final.to_excel(DF_OUT,index=None)
