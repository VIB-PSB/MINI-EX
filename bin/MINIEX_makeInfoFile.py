
import pandas,collections,seaborn, sys


MAT=sys.argv[1]
GRNB_OUT=sys.argv[2]
MOTENR_OUT=sys.argv[3]
TFFILT_OUT=sys.argv[4]
INFO_TF=sys.argv[5]
CELLS=sys.argv[6]
IDS=sys.argv[7]
OUT=sys.argv[8]


matrix=pandas.read_csv(MAT,sep='\t',index_col=0)

exp_genes = list(matrix.index)

cell2type={}
with open(CELLS) as f:
    for line in f:
        cell2type[line.strip().split('\t')[0]]=str(line.strip().split('\t')[1])
cell2type_counter=collections.Counter(list(cell2type.values()))  

tfs_all=list(pandas.read_csv(INFO_TF, header=None, usecols=[0], sep='\t')[0])

matrix=matrix[matrix.index.isin(tfs_all)]  #subset matrix to only keep TFs
matrix=matrix.transpose()
matrix=matrix[matrix.index.isin(cell2type.keys())] #subset matrix to only keep annotated cells

tfs_exp=list(set(exp_genes)& set(tfs_all))

tfs_grnb=list(set(pandas.read_csv(GRNB_OUT, header=None, usecols=[0], sep='\t')[0]))

tfs_enrich=list(set(pandas.read_csv(MOTENR_OUT, header=None, usecols=[0], sep='\t')[0]))

tfs_filt=list(set(pandas.read_csv(TFFILT_OUT, header=None, usecols=[0], sep='\t')[0]))

df_info=pandas.DataFrame(0,index=sorted(tfs_all),columns=['isTF_expressed','isTF_expressionRegulons','isTF_TFBSRegulons','isTF_finalRegulons'])
for tf in tfs_all:
    if tf in tfs_exp:
        df_info['isTF_expressed'][tf]=1
    if tf in tfs_grnb:
        df_info['isTF_expressionRegulons'][tf]=1
    if tf in tfs_enrich:
        df_info['isTF_TFBSRegulons'][tf]=1
    if tf in tfs_filt:
        df_info['isTF_finalRegulons'][tf]=1

identities={}
with open(IDS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        identities[str(spl[0])]=spl[1]
        
columns=[identities[i]+'-'+i for i in cell2type_counter.keys()]        
df_perc=pandas.DataFrame(0.0,index=tfs_exp,columns=sorted(columns))

    
for tf in tfs_exp:
    cell_tot = collections.Counter([cell2type[i] for i in matrix[matrix[tf]>1].index.tolist() ]) #in how many cells the TF is exp
    temp=[]
    for ct in cell_tot:
        percent=round((cell_tot[ct]/cell2type_counter[ct])*100,1)
        df_perc[identities[ct]+'-'+ct][tf]=percent

        
df_all=df_info.join(df_perc)
df_all.to_csv(OUT,sep='\t')
