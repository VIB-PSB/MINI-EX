"""
Identifies whether to use the standard or the reference Borda ranking.

From the provided list of regulons, extracts their TFs and identifies the list
of relevant TFs (i.e., having at least one entry in the GO annotation file
containing in its description at least one term of interest).

If at least two relevant TFs are found: use the reference procedure.
Otherwise: use the standard procedure (which will also take into account
the GO enrichment value).
"""

import pandas as pd
import sys
             
REGULONS_FILE          = sys.argv[1]
GO_ANNOTATIONS_FILE    = sys.argv[2]
TERMS_OF_INTEREST_FILE = sys.argv[3]


if TERMS_OF_INTEREST_FILE.startswith('.dummy_path'):
    # if terms of interest aren't provided: use the standard procedure
    print("std")
else:
    # collect TFs from the list of regulons (first column)
    all_tfs = set(pd.read_csv(REGULONS_FILE, sep='\t', header=None)[0].unique())
    
    terms_of_interest = pd.read_csv(TERMS_OF_INTEREST_FILE, sep='\t', header=None)[0].unique().tolist()

    # collect list of relevant TFs: are present in the GO annotations file 
    # and their description contains at least one term of interest

    go_annotations_df = pd.read_csv(GO_ANNOTATIONS_FILE, sep='\t', header=None, names=['go_term', 'gene_id', 'evidence_code', 'go_description'])

    # a term is relevant if either its term or description contain at least one term of interest
    pattern = '|'.join(terms_of_interest)
    go_annotations_df['term_is_relevant'] = go_annotations_df[['go_term', 'go_description']].apply(lambda x: x.str.contains(pattern, case=False, na=False)).any(axis=1)
    releavant_annotated_genes = set(go_annotations_df[go_annotations_df['term_is_relevant']]['gene_id'])
    relevant_tfs = all_tfs.intersection(releavant_annotated_genes)

    # if at least two relevant TFs were found: use "ref" approach, otherwise: use "std" procedure
    if len(relevant_tfs) < 2:
        print("std")
    else:
        print("ref")
