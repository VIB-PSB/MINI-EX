import sys
             
REGULONS=sys.argv[1]
GOS=sys.argv[2]

OUT_SET=sys.argv[3]
           
OUT_FEAT=sys.argv[4]

reg2genes={}
allgenes=[]
with open(REGULONS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        tf=spl[0]
        allgenes.append(tf)
        allgenes+=spl[2].rsplit(',')
        if tf+'_'+spl[1] in reg2genes:
            reg2genes[tf+'_'+spl[1]]+=spl[2].rsplit(',')
            reg2genes[tf+'_'+spl[1]]+=[tf]
        else:
            reg2genes[tf+'_'+spl[1]]=spl[2].rsplit(',')
            reg2genes[tf+'_'+spl[1]]+=[tf]

allgenes=list(set(allgenes))

foutS=open(OUT_SET,'w')
for r in reg2genes:
    reg2genes[r]=list(set(reg2genes[r]))
    for g in reg2genes[r]:
        foutS.write(r+'\t'+g+'\n')
        
foutS.close()


foutF=open(OUT_FEAT,'w')        
with open(GOS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        if spl[1] in allgenes:
            foutF.write(spl[0]+'\t'+spl[1]+'\n')
            
foutF.close()