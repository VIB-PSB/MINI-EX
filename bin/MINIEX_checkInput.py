import sys
import pandas as pd

EXPRESSION_MATRIX=sys.argv[1].split(' ')
MARKERS_OUT=sys.argv[2].split(' ')
CELLS2CLUSTERS=sys.argv[3].split(' ')
CLUSTER2IDENT=sys.argv[4].split(' ')
TF_LIST=sys.argv[5].split(' ')
TERMS_OF_INTEREST=sys.argv[6].split(' ')
GRNBOOST_OUT=sys.argv[7].split(' ')
FEATURE_FILE_MOTIFS=sys.argv[8].split(' ')
INFO_TF=sys.argv[9].split(' ')
GO_FILE=sys.argv[10].split(' ')
ALIAS=sys.argv[11].split(' ')
DO_MOTIF_ANALYSIS=sys.argv[12]
TOP_MARKERS=sys.argv[13]
EXPRESSION_FILTER=sys.argv[14]
MOTIF_FILTER=sys.argv[15]
TOP_REGULONS=sys.argv[16]

########## NUMERICAL PARAMETERS ##########
def value_is_a_positive_integer(value):
    try:
        int_value = int(value)
        return int_value >= 0
    except ValueError:
        return False

if not value_is_a_positive_integer(TOP_MARKERS):
    raise Exception(f"Incorrect value provided for the parameter 'topMarkers': '{TOP_MARKERS}'! It can only have positive integer values.")

if not value_is_a_positive_integer(TOP_REGULONS):
    raise Exception(f"Incorrect value provided for the parameter 'topRegulons': '{TOP_REGULONS}'! It can only have positive integer values.")

if DO_MOTIF_ANALYSIS != "true" and DO_MOTIF_ANALYSIS != "false":
    raise Exception(f"Incorrect value provided for the parameter 'doMotifAnalysis': '{DO_MOTIF_ANALYSIS}'. Only boolean values are allowed.")

if MOTIF_FILTER != "TF-F_motifs" and MOTIF_FILTER != "TF_motifs":
    raise Exception(f"Incorrect value provided for the parameter 'motifFilter': '{MOTIF_FILTER}'! The allowed values are 'TF-F_motifs' or 'TF_motifs'.")

if not value_is_a_positive_integer(EXPRESSION_FILTER) or int(EXPRESSION_FILTER) > 100:
    raise Exception(f"Incorrect value provided for the parameter 'expressionFilter': '{EXPRESSION_FILTER}'! The allowed values are between 0 and 100.")


########## PATH PARAMETERS ##########

def path_is_dummy(path): # checks whether the user specified "null" as path to that file
    return path[0].startswith('dummy_path')

def count_lines(file_path):  # counts number of lines in the provided file without reading it
    with open(file_path, 'r') as file:
        line_count = sum(1 for line in file)
    return line_count

# collecting names of all input files
all_dataset_file_names = EXPRESSION_MATRIX + MARKERS_OUT + CELLS2CLUSTERS + CLUSTER2IDENT
all_file_names = all_dataset_file_names + TF_LIST + INFO_TF + ALIAS
if not path_is_dummy(TERMS_OF_INTEREST):
    all_file_names += TERMS_OF_INTEREST
if not path_is_dummy(GRNBOOST_OUT):
    all_file_names += GRNBOOST_OUT
if not path_is_dummy(FEATURE_FILE_MOTIFS):
    all_file_names += FEATURE_FILE_MOTIFS
if not path_is_dummy(GO_FILE):
    all_file_names += GO_FILE


#### GENERAL VERIFICATIONS ####

# CHECK: if the GO file is null, then terms of interest is also null
if path_is_dummy(GO_FILE):
    if not path_is_dummy(TERMS_OF_INTEREST):
        raise Exception("When the GO file is set to 'null', the terms of interest file should also be set to null!")
        
# CHECK: equal number of files was provided for each file category, one per dataset
if len(EXPRESSION_MATRIX) != len(MARKERS_OUT):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of markers files ({len(MARKERS_OUT)})!")
if len(EXPRESSION_MATRIX) != len(CELLS2CLUSTERS):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of cell2clusters files ({len(CELLS2CLUSTERS)})!")
if len(EXPRESSION_MATRIX) != len(CLUSTER2IDENT):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of identities files ({len(CLUSTER2IDENT)})!")
if not path_is_dummy(GRNBOOST_OUT) and len(EXPRESSION_MATRIX) != len(GRNBOOST_OUT):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of GRNBoost output files ({len(GRNBOOST_OUT)})!")


#### DATASET-RELATED VERIFICATIONS ####

# as MINI-EX can handle multiple datasets at the same time, the list of datasets is extracted
# and the following verifications are performed for each dataset individually

data_sets = [file_name.split('_')[0] for file_name in EXPRESSION_MATRIX]

# will collect statistics about the input files for each dataset
stats_df = pd.DataFrame(columns=['dataset', 'cells', 'genes', 'clusters', 'tissues'])
stats_df['dataset'] = data_sets
stats_df.set_index('dataset', inplace=True)

# CHECK: all dataset files should start with the name of dataset followed by an underscore
if (not all('_' in file_name for file_name in all_dataset_file_names) or
    (GRNBOOST_OUT != None and not all('_' in file_name for file_name in GRNBOOST_OUT))):
    raise Exception("All input files must contain an underscore!")

# performs checks common to all input dataset files
def check_dataset_file(file_name: str):
    # CHECK: no quotes allowed in input files
    with open(file_name) as a_file:
        if '"' in a_file.read():
            raise Exception(f"Quotes (\") detected in file '{file_name}'!")
    try:
        df = pd.read_csv(file_name, delimiter='\t', header=None, index_col=0)
        row_count = len(df)
    except pd.errors.ParserError: # in some files generated from Seurat, the header line has one column missing -> skip it
        df = pd.read_csv(file_name, delimiter='\t', header=None, index_col=0, skiprows=1)
        row_count = len(df) + 1
    finally:
        # CHECK: no trailing newlines present in the input files
        if row_count != count_lines(file_name):
            raise Exception(f"Empty lines found in '{file_name}'! ({len(df)} vs {count_lines(file_name)})")
        # CHECK: values in the first column of each data file must be unique
        # the only exception is the Seurat markers file, so it is not checked
        if df.index.has_duplicates and file_name not in MARKERS_OUT:
            raise Exception(f"Duplicated indexes detected in file '{file_name}'!")
    
for file_name in all_dataset_file_names:
    check_dataset_file(file_name)

for data_set in data_sets:
    # CHECK: each dataset must have the four input files (or five, if GRNBoost output files are provided)
    prefix = f"{data_set}_"
    if not any(file_name.startswith(prefix) for file_name in MARKERS_OUT):
        raise Exception(f"No markers file provided for the dataset {data_set}!")
    if not any(file_name.startswith(prefix) for file_name in CELLS2CLUSTERS):
        raise Exception(f"No cells2clusters file provided for the dataset {data_set}!")
    if not any(file_name.startswith(prefix) for file_name in CLUSTER2IDENT):
        raise Exception(f"No identities file provided for the dataset {data_set}!")
    if not path_is_dummy(GRNBOOST_OUT) and not any(file_name.startswith(prefix) for file_name in GRNBOOST_OUT):
        raise Exception(f"No GRNBoost output file provided for the dataset {data_set}!")
    
    # read the four input files for the current dataset as pandas dataframes
    cells2clusters_file_name = [file_name for file_name in CELLS2CLUSTERS if file_name.startswith(prefix)][0]
    identities_file_name = [file_name for file_name in CLUSTER2IDENT if file_name.startswith(prefix)][0]
    matrix_file_name = [file_name for file_name in EXPRESSION_MATRIX if file_name.startswith(prefix)][0]
    markers_file_name = [file_name for file_name in MARKERS_OUT if file_name.startswith(prefix)][0]
    cells2clusters_df = pd.read_csv(cells2clusters_file_name,delimiter='\t',names=['cell_id','cluster_id','cluster_annotation'],dtype={'cell_id':str,'cluster_id':str,'cluster_annotation':str})[['cluster_id']].drop_duplicates().sort_values(by='cluster_id').reset_index(drop=True)
    identities_df = pd.read_csv(identities_file_name,delimiter='\t',usecols=[0,1],names=['cluster_id','cluster_annotation'],dtype={'cluster_id':str,'cluster_annotation':str}).drop_duplicates().sort_values(by='cluster_id').reset_index(drop=True)
    matrix_df = pd.read_csv(matrix_file_name,delimiter='\t',header=0,index_col=0)
    markers_df = pd.read_csv(markers_file_name,delimiter='\t',skiprows=1,names=['gene_id', 'p_val', 'avg_logFC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster', 'gene'])

    # CHECK: no underscores allowed in cluster names
    if any('_' in cluster_name for cluster_name in identities_df['cluster_annotation'].tolist()):
        raise Exception(f"Underscores detected in cluster names of the '{data_set}' data set!")
    
    # CHECK: cluster id and identities must correspond between cells2clusters and cluster2ident files
    identities_cluster_id_df = identities_df[['cluster_id']].drop_duplicates().sort_values(by='cluster_id').reset_index(drop=True).copy()
    if not cells2clusters_df.equals(identities_cluster_id_df):
        raise Exception(f"Cluster identities differ between cells2cluster and identities files for the '{data_set}' data set!")

    # CHECK: only up-regulated markers should be present
    if (markers_df['avg_logFC'] < 0).any():
        raise Exception(f"Down-regulated markers detected for the '{data_set}' data set!")

    # collect statistics
    stats_df.loc[data_set, 'cells'] = matrix_df.shape[1]
    stats_df.loc[data_set, 'genes'] = matrix_df.shape[0]
    stats_df.loc[data_set, 'clusters'] = len(identities_df.index)
    stats_df.loc[data_set, 'tissues'] = len(identities_df[['cluster_annotation']].drop_duplicates().index)


# if the process is here, that means that all the tests passed, otherwise an exception is raized and the process is interrupted
print("")
print("== INPUT VALIDATION ==========================================")
print("Input files passed validation tests. Retrieved following data:")
print(stats_df)
print("")
print("== INPUT FILES ===============================================")
print(f"Expression matrix file(s) : {' / '.join(EXPRESSION_MATRIX)}")
print(f"Seurat markers file(s)    : {' / '.join(MARKERS_OUT)}")
print(f"Cells to clusters file(s) : {' / '.join(CELLS2CLUSTERS)}")
print(f"Cluster identities file(s): {' / '.join(CLUSTER2IDENT)}")
print(f"GRNBoost output file(s)   : {' / '.join(GRNBOOST_OUT) if not path_is_dummy(GRNBOOST_OUT) else 'NOT PROVIDED'}")
print(f"Transcription factor file : {' / '.join(TF_LIST)}")
print(f"TF info file              : {' / '.join(INFO_TF)}")
print(f"Gene aliases file         : {' / '.join(ALIAS)}")
print(f"Motifs feature file       : {' / '.join(FEATURE_FILE_MOTIFS) if not path_is_dummy(FEATURE_FILE_MOTIFS) else 'NOT PROVIDED'}")
print(f"GO file                   : {' / '.join(GO_FILE) if not path_is_dummy(GO_FILE) else 'NOT PROVIDED'}")
print(f"Terms of interest file    : {' / '.join(TERMS_OF_INTEREST) if not path_is_dummy(TERMS_OF_INTEREST) else 'NOT PROVIDED'}")
print("")
print("== MINI-EX PARAMETERS ========================================")
if not path_is_dummy(TERMS_OF_INTEREST):
    terms_of_interest_df = pd.read_csv(TERMS_OF_INTEREST[0], names=['term'])
    print(f"Terms of interest         : {' / '.join(terms_of_interest_df['term'])}")
print(f"doMotifAnalysis           : {DO_MOTIF_ANALYSIS}")
print(f"topMarkers                : {TOP_MARKERS}")
print(f"expressionFilter          : {EXPRESSION_FILTER}")
print(f"motifFilter               : {MOTIF_FILTER}")
print(f"topRegulons               : {TOP_REGULONS}")
print("")