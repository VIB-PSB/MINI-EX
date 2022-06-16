
import pandas,collections,seaborn, sys

# MAIN='//psb.ugent.be/shares/research/deepseq/ngsproject_scgrn/results/cafer/miniex_osa/'
MAT=sys.argv[1]#MAIN+'files/OsaGRN_matrix.txt'
GRNB_OUT=sys.argv[2]#MAIN+'files/OsaGRN_grnboost2.txt'
MOTENR_OUT=sys.argv[3]#MAIN+'work/5e/14d5fdb00c9d6fc66daae0fea8ad1c/OsaGRN_enrichedRegulons.txt'
TFFILT_OUT=sys.argv[4]#MAIN+'regulons_output/OsaGRN_regulons.txt'
INFO_TF=sys.argv[5]#'//psb.ugent.be/shares/biocomp/groups/group_cig/cafer/vib-psb-github/MINI-EX/data_osa/osa_TF2fam2mot.txt'#sys.argv[2]
CELLS=sys.argv[6]#MAIN+'files/OsaGRN_cells2clusters.txt'
IDS=sys.argv[7]#MAIN+'files/OsaGRN_identities.txt'
OUT=sys.argv[8]


matrix=pandas.read_csv(MAT,sep='\t',index_col=0)

exp_genes = list(matrix.index)
tfs_all=list(pandas.read_csv(INFO_TF, header=None, usecols=[0], sep='\t')[0])

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
        
        
        
cell2type={}
with open(CELLS) as f:
    for line in f:
        cell2type[line.strip().split('\t')[0]]=line.strip().split('\t')[1]

cell2type_counter=collections.Counter(list(cell2type.values()))        
   

identities={}
with open(IDS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        identities[spl[0]]=spl[1]
        
columns=[identities[i]+'-'+i for i in cell2type_counter.keys()]        
df_perc=pandas.DataFrame(0.0,index=tfs_exp,columns=sorted(columns))

    
matrix=matrix.transpose()

for tf in tfs_exp:
    cell_tot = collections.Counter([cell2type[i] for i in matrix[matrix[tf]>1].index.tolist() ]) #in how many cells the TF is exp
    temp=[]
    for ct in cell_tot:
        percent=round((cell_tot[ct]/cell2type_counter[ct])*100,1)
        df_perc[identities[ct]+'-'+ct][tf]=percent

        
df_all=df_info.join(df_perc)
df_all.to_csv(OUT,sep='\t')
