
"""
Identifies whether to use the standard Borda ranking or also take into account
GO enrichment q-values.

From the provided list of regulons, extracts their TFs and identifies the list
of relevant TFs (i.e., having at least one entry in the GO annotation file
containing in its description at least one term of interest).

If at least two relevant TFs are found: GO enrichment q-value should be taken into account.
Otherwise: standard Borda ranking should be computed.
"""

import sys,pandas
             
REGULONS=sys.argv[1]  # path to list of regulons
GOS=sys.argv[2]  # path to GO annotations
GOI=sys.argv[3]  # path to list of terms of interest

# collect TFs from the list of regulons
allTFs=[]
with open(REGULONS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        allTFs+=[spl[0]]
allTFs=list(set(allTFs))

# list of terms of interest
listOfterms=list(set(pandas.read_csv(GOI,sep='\t',header=None)[0]))

# collect list of relevant TFs: are present in GO annotations file 
# and their description contains at least one term of interest
myTFs=[]    
with open(GOS) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')  # columns are: GO term - gene - evidence code - annotation
        if any(s in spl[3] for s in listOfterms) and spl[1] in allTFs:
            myTFs.append(spl[1])      
myTFs=list(set(myTFs))            

# if at least two relevant TFs were found: use "ref" approach, otherwise: use "std" approach
if len(myTFs)<2:
    print('std')
else:
    print('ref')
