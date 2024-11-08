
import pandas, collections, sys

###define which cell-types the cluster belong to
CELLID=sys.argv[1]      
REGULONS=sys.argv[2]
FILE_ALIAS=sys.argv[3]
TOPS=sys.argv[4]
GOS=sys.argv[5]
GOI=sys.argv[6]
CENTRALITIES=sys.argv[7]
ENRICHMENT=sys.argv[8]
ALLMARKERS=sys.argv[9]
OUTFILE=sys.argv[10]

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
        geneID, geneSymbol = line.rstrip().split('\t')[0], line.rstrip().split('\t')[1]
        if geneID in allTFs:
            if geneID not in aliases:
                aliases[geneID] = [geneSymbol]
            else:
                aliases[geneID].append(geneSymbol)

for geneID, symbol in aliases.items():
    aliases[geneID] = "/ ".join(aliases[geneID]) # merge the final list of aliases (gene symbols) into one string

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
termsOI,termsNOI={},{}
genesOI,genesNOI=[],[]
gene2go={}

listOfterms=list(set(pandas.read_csv(GOI,sep='\t',header=None)[0]))

with open(GOS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        if any(s in spl[3] for s in listOfterms):
            termsOI[spl[0]]=spl[3]
            genesOI.append(spl[1])
        else:
            termsNOI[spl[0]]=spl[3]
            genesNOI.append(spl[1])
        if spl[1] in gene2go:
            gene2go[spl[1]]+=[spl[0]]
        else:
            gene2go[spl[1]]=[spl[0]]

genesOI=list(set(genesOI)) 
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


###for root GO enrichment (among the TGs) take the term with lowest enrichment (below)
dic_go_tmp={}
with open(ENRICHMENT) as f:
    for line in f:
        if line.startswith('#') or line.startswith('set_id'): # skip the comments and the header
            pass
        else:
            spl=line.rstrip().rsplit('\t')
            if spl[1] in termsOI:
                if spl[0] in dic_go_tmp:
                    dic_go_tmp[spl[0]]+=[[float(spl[3]),spl[1],termsOI[spl[1]],int(spl[8])]]
                else:
                    dic_go_tmp[spl[0]]=[[float(spl[3]),spl[1],termsOI[spl[1]],int(spl[8])]]

dic_go={}
for el in dic_go_tmp:
    if len(dic_go_tmp[el])==1:
        dic_go[el]=dic_go_tmp[el][0]
    else:
        tmp=[]
        for me in dic_go_tmp[el]:
            tmp.append(me[0])
        mino=min(tmp)
        dic_go[el]=[m for m in dic_go_tmp[el] if m[0]==mino][0]
        
###add combos (gene_cluster) which are not enriched and give nan vals                    
for el in dic_enirch:
    if el in dic_go:
        pass
    else:
        dic_go[el]=[float('nan'),'-','-',float('nan')]

# create a dictionary of clusters vs their DEGs
markers_df = pandas.read_csv(ALLMARKERS, sep='\t',dtype={'cluster': str})
markers_df = markers_df[markers_df['cluster'].isin(cellTyp_mtx.keys())]
markers_df['cluster'] = markers_df['cluster'].apply(lambda x: f"Cluster_{x}")
markers_df = markers_df[markers_df['p_val_adj'] <= 0.05]
clu2gene = markers_df.groupby('cluster')['gene'].unique().apply(list).to_dict()
    
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
    if geneid in genesOI:
        whichGO=sorted([goo for goo in gene2go[geneid] if goo in termsOI])
        whichGOdesc=[termsOI[goo] for goo in whichGO]
        df.append([geneid,aliases[geneid],'relevant_known_TF',','.join(whichGO),','.join(whichGOdesc),cellTyp_mtx[clusterid.replace('Cluster_','')]+'_'+clusterid,cellTyp_mtx[clusterid.replace('Cluster_','')],tfclu2de[ele],regPerClu[clusterid],dic_reg[ele],dic_enirch[ele]]+dic_centrality[ele]+dic_go[ele])
    elif geneid in genesNOI:
        whichGO=sorted([goo for goo in gene2go[geneid] if goo in termsNOI])
        whichGOdesc=[termsNOI[goo] for goo in whichGO]
        df.append([geneid,aliases[geneid],'known_TF',','.join(whichGO),','.join(whichGOdesc),cellTyp_mtx[clusterid.replace('Cluster_','')]+'_'+clusterid,cellTyp_mtx[clusterid.replace('Cluster_','')],tfclu2de[ele],regPerClu[clusterid],dic_reg[ele],dic_enirch[ele]]+dic_centrality[ele]+dic_go[ele])
    else:
        df.append([geneid,aliases[geneid],'unknown_TF','-','-',cellTyp_mtx[clusterid.replace('Cluster_','')]+'_'+clusterid,cellTyp_mtx[clusterid.replace('Cluster_','')],tfclu2de[ele],regPerClu[clusterid],dic_reg[ele],dic_enirch[ele]]+dic_centrality[ele]+dic_go[ele])

df=pandas.DataFrame(df)   
df.columns=['TF','alias','hasTFrelevantGOterm','GOterm','GOdescription','cluster','celltype','isTF_DE','totRegInCluster','#TGs','qval_cluster','out-degree','closeness','betweenness','GO_enrich_qval','GO_enrich_term','GO_enrich_desc','#TGs_withGO']

df.to_csv(OUTFILE,sep='\t',index=None)
