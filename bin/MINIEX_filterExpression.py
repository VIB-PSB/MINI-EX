import sys,pandas, collections

MATRIX=sys.argv[1]
TF_FILE=sys.argv[2]

CELLS=sys.argv[3]

     
matrix=pandas.read_csv(MATRIX,sep='\t',index_col=0)
cells_inMat=list(matrix.columns)
tf_list=pandas.read_csv(TF_FILE,sep='\t',header=None)[0].to_list()

matrix=matrix[matrix.index.isin(tf_list)] #subset matrix to only keep TFs


cell2type={}
with open(CELLS) as f:
    for line in f:
        if line.strip().split('\t')[0] in cells_inMat:
            cell2type[line.strip().split('\t')[0]]=str(line.strip().split('\t')[1])

cell2type_counter=collections.Counter(list(cell2type.values()))        

matrix=matrix.transpose()

matrix=matrix[matrix.index.isin(cell2type.keys())] #subset matrix to only keep annotated cells
matrix_col=list(matrix.columns)

tfsInMatrix=list(set(matrix_col)&set(tf_list))


   
perc=sys.argv[4]   
tf2clu={}

for tf in tfsInMatrix:
    cell_tot = collections.Counter([cell2type[i] for i in matrix[matrix[tf]>1].index.tolist() ]) #in how many cells the TF is exp
    temp=[]
    for ct in cell_tot:
        if cell_tot[ct]>= round((int(perc)/100)*cell2type_counter[ct]):
            temp.append(ct)
    if len(temp) != 0:
        tf2clu[tf]=temp
    else:
        tf2clu[tf]='none'



def filtering(inf,out):      
    fin=open(inf,'r')
    fout=open(out,'w')
    dic={}
    for line in fin:
        if line.startswith('#') or line.startswith('set_id'): # skip the comments and the header
            pass
        else:
            spl=line.rstrip().rsplit('\t')
            reg=spl[0]
            target_genes=spl[9].rsplit(",")
            ct=spl[1].replace('Cluster_','')
            if ct in  tf2clu[reg]:
                dic[reg+'_'+ct]=target_genes
                
    for k in dic:
        dic[k]=list(set(dic[k]))
                        
        splK=k.rstrip().rsplit('_')
        
        fout.write('_'.join(splK[:-1])+'\tCluster_'+splK[-1]+'\t'+','.join(dic[k])+'\n')
    fout.close()
        
infi = sys.argv[5]
outi = sys.argv[6]
filtering(infi,outi)