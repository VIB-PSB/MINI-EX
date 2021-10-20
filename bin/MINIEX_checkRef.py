
import sys,pandas
             
REGULONS=sys.argv[1]        
GOS=sys.argv[2]
GOI=sys.argv[3]

allTFs=[]
with open(REGULONS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        allTFs+=[spl[0]]


allTFs=list(set(allTFs))

listOfterms=list(set(pandas.read_csv(GOI,sep='\t',header=None)[0]))

myTFs=[]    
with open(GOS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        if any(s in spl[3] for s in listOfterms) and spl[1] in allTFs:
            myTFs.append(spl[1])

            
myTFs=list(set(myTFs))            

if len(myTFs)<2:
    print('std')
else:
    print('ref')
