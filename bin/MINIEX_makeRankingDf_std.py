import pandas, collections, sys

###define which cell-types the cluster belong to
CELLID=sys.argv[1]  
REGULONS=sys.argv[2]
FILE_ALIAS=sys.argv[3]
TOPS=sys.argv[4]
GOS=sys.argv[5]
CENTRALITIES=sys.argv[6]
ALLMARKERS=sys.argv[7]
OUTFILE=sys.argv[8]

cellTyp_mtx = {}
for line in open(CELLID):
    cellTyp_mtx[str(line.rstrip().rsplit('\t')[0])]=line.rstrip().rsplit('\t')[1]
                                           
###load final regulons (filtered for TF being expressed in at least 5% of cells in the cluster), #TGs and #TGswithGWAS
dic_reg={}
with open(REGULONS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        dic_reg[spl[0]+'_'+spl[1]]=len(spl[2].rsplit(','))
allTFs=list(set(['_'.join(i.rsplit('_')[:-2]) for i in list(dic_reg.keys())]))   

###define aliases for TFs
aliases={}
with open(FILE_ALIAS) as f:
    for line in f:
        if line.rstrip().rsplit('\t')[0] in allTFs:
            aliases[line.rstrip().rsplit('\t')[0]]=line.rstrip().rsplit('\t')[1]

for tf in allTFs:
    if tf not in aliases:
        aliases[tf]=tf
regPerClu=collections.Counter(['_'.join(i.rsplit('_')[-2:]) for i in dic_reg])

###load output regulons to get qvals corresponding to regulons enrichment in clusters
dic_enirch={}
with open(TOPS) as f:
    for line in f:
        if line.startswith('#'):
            pass
        else:
            spl=line.rstrip().rsplit('\t')
            if spl[0]+'_'+spl[1] in dic_reg:
                dic_enirch[spl[0]+'_'+spl[1]]=float(spl[3])
            
###do the same s above for non root GOs
termsNOI={}
genesNOI=[]
gene2go={}

with open(GOS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        termsNOI[spl[0]]=spl[3]
        genesNOI.append(spl[1])
        if spl[1] in gene2go:
            gene2go[spl[1]]+=[spl[0]]
        else:
            gene2go[spl[1]]=[spl[0]]

genesNOI=list(set(genesNOI))


###load centrality measurments for regulons
dic_centrality={}
with open(CENTRALITIES) as f:
    header=f.readlines()[0].rstrip().rsplit('\t')
with open(CENTRALITIES) as f:
    for line in f.readlines()[1:]:
        spl=line.rstrip().rsplit('\t')
        for j in range(len(header)):
            clu='_'.join(header[j].rsplit('_')[1:3])
            if spl[0]+'_'+clu in dic_reg:
                wher=[i for i in range(len(header)) if header[i].endswith(clu)]
                dic_centrality[spl[0]+'_'+clu]=[float(spl[wher[0]]),float(spl[wher[1]]),float(spl[wher[2]])]
            else: ###add zeros for the gene_cluster for which measurements are not done
                if spl[0]+'_'+clu in dic_enirch:
                    dic_centrality[spl[0]+'_'+clu]=[0.0,0.0,0.0]



clusters=["Cluster_"+i for i in cellTyp_mtx.keys()]
clu2gene={}	
for line in open(ALLMARKERS,'r').readlines()[1:]:
    spl=line.rstrip().rsplit('\t')
    if float(spl[5]) <=0.05:# and float(spl[2])>=0.5 and float(spl[3]) >= 0.5:
        if 'Cluster_'+spl[6] in clusters:
            if 'Cluster_'+spl[6] in clu2gene:
                clu2gene['Cluster_'+spl[6]]+=[spl[7]]
            else:
                clu2gene['Cluster_'+spl[6]]=[spl[7]]

for clu in clu2gene:
    clu2gene[clu]=list(set(clu2gene[clu]))
    
tfclu2de={}
for el in dic_reg:
    if '_'.join(el.rsplit('_')[:-2]) in clu2gene["Cluster_"+el.rsplit('_')[-1]]:
        tfclu2de[el]=1
    else:
        tfclu2de[el]=0

df=[]
for ele in dic_enirch:
    geneid='_'.join(ele.rsplit('_')[:-2])
    clusterid='_'.join(ele.rsplit('_')[-2:])
    if geneid in genesNOI:
        whichGO=sorted([goo for goo in gene2go[geneid] if goo in termsNOI])
        whichGOdesc=[termsNOI[goo] for goo in whichGO]
        df.append([geneid,aliases[geneid],'known_TF',','.join(whichGO),','.join(whichGOdesc),cellTyp_mtx[clusterid.replace('Cluster_','')]+'_'+clusterid,cellTyp_mtx[clusterid.replace('Cluster_','')],tfclu2de[ele],regPerClu[clusterid],dic_reg[ele],dic_enirch[ele]]+dic_centrality[ele])
    else:
        df.append([geneid,aliases[geneid],'unknown_TF','-','-',cellTyp_mtx[clusterid.replace('Cluster_','')]+'_'+clusterid,cellTyp_mtx[clusterid.replace('Cluster_','')],tfclu2de[ele],regPerClu[clusterid],dic_reg[ele],dic_enirch[ele]]+dic_centrality[ele])

df=pandas.DataFrame(df)   
df.columns=['TF','alias','hasTFrelevantGOterm','GOterm','GOdescription','cluster','celltype','isTF_DE','totRegInCluster','#TGs','qval_cluster','out-degree','closeness','betweenness']

df.to_csv(OUTFILE,sep='\t',index=None)
