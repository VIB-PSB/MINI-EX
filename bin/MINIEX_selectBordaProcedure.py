"""
Identifies whether to use the standard or the reference Borda ranking.

From the provided list of regulons, extracts their TFs and identifies the list
of relevant TFs (i.e., having at least one entry in the GO annotation file
containing in its description at least one term of interest).

If at least two relevant TFs are found: use the reference procedure.
Otherwise: use the standard procedure (which will also take into account
the GO enrichment value).
"""

import sys,pandas
             
REGULONS_FILE=sys.argv[1]
GO_ANNOTATIONS_FILE=sys.argv[2]
TERMS_OF_INTEREST_FILE=sys.argv[3]

# collect TFs from the list of regulons
all_tfs = []
with open(REGULONS_FILE) as f:
    for line in f:
        spl=line.rstrip().rsplit('\t')
        all_tfs += [spl[0]]
all_tfs = list(set(all_tfs))


if TERMS_OF_INTEREST_FILE.startswith('.dummy_path'):
    # if terms of interest aren't provided: use the standard procedure
    print("std")
else:
    terms_of_interest = list(set(pandas.read_csv(TERMS_OF_INTEREST_FILE,sep='\t',header=None)[0]))

    # collect list of relevant TFs: are present in the GO annotations file 
    # and their description contains at least one term of interest
    relevant_tfs = []
    with open(GO_ANNOTATIONS_FILE) as f:
        for line in f:
            spl=line.rstrip().rsplit('\t')  # columns are: GO term - gene - evidence code - GO description
            gene_id = spl[1]
            go_description = spl[3]
            if any(term_of_interest in go_description for term_of_interest in terms_of_interest) and gene_id in all_tfs:
                relevant_tfs.append(gene_id)
    myTFs=list(set(relevant_tfs))   

    # if at least two relevant TFs were found: use "ref" approach, otherwise: use "std" procedure
    if len(myTFs)<2:
        print('std')
    else:
        print('ref')
