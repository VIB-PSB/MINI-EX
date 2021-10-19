import pandas, sys
import networkx as nx
       
REGULONS=sys.argv[1]
OUTFILE=sys.argv[2]
dic={}
clust=[]
with open(REGULONS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        clust.append(spl[1])
        if spl[0]+'_'+spl[1] in dic:
            dic[spl[0]+'_'+spl[1]]+=spl[2].rsplit(',')
        else:
            dic[spl[0]+'_'+spl[1]]=spl[2].rsplit(',')
          
clust=list(set(clust))    
        
doi={}
for c in clust:
    doi2={}
    for t in dic:
        if c ==  '_'.join(t.rsplit('_')[-2:]) and '_'.join(t.rsplit('_')[:-2]) in doi2:
            doi2['_'.join(t.rsplit('_')[:-2])]+=list(set(dic[t]))
        elif c ==  '_'.join(t.rsplit('_')[-2:]) and '_'.join(t.rsplit('_')[:-2]) not in doi2:
            doi2['_'.join(t.rsplit('_')[:-2])]=list(set(dic[t]))
        doi[c]=doi2

cols=[]
for c in sorted(doi.keys()):
    cols+='degout_'+c,'clos_'+c,'bet_'+c
indi=list(set(['_'.join(i.rsplit('_')[:-2]) for i in dic.keys()]))
df=pandas.DataFrame(float('nan'),index=indi,columns=cols)

for c in doi:
    print(c)
    G = nx.DiGraph()
    allEdges=[]
    for ttt in doi[c]: 
        listOfEdges=[]
        for gg in doi[c][ttt]:
            allEdges.append((ttt,gg))
            listOfEdges.append((ttt,gg))

    G.add_edges_from(allEdges)
    deg_out=nx.out_degree_centrality(G)
    clos=nx.closeness_centrality(G)
    bet=nx.betweenness_centrality(G)

    for tt in doi[c]:
        df['degout_'+c][tt]=deg_out[tt]
        df['clos_'+c][tt]=clos[tt]
        df['bet_'+c][tt]=bet[tt]
      
df.to_csv(OUTFILE,sep='\t')
