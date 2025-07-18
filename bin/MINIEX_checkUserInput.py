import subprocess
import sys
import pandas as pd

EXPRESSION_MATRIX=sys.argv[1].split(' ')
MARKERS_OUT=sys.argv[2].split(' ')
CELLS2CLUSTERS=sys.argv[3].split(' ')
CLUSTER2IDENT=sys.argv[4].split(' ')
TF_LIST=sys.argv[5]
TERMS_OF_INTEREST=sys.argv[6]
GRNBOOST_OUT=sys.argv[7].split(' ')
FEATURE_FILE_MOTIFS=sys.argv[8]
INFO_TF=sys.argv[9]
GO_FILE=sys.argv[10]
ALIAS=sys.argv[11]
ENRICHMENT_BACKGROUND=sys.argv[12]
DO_MOTIF_ANALYSIS=sys.argv[13]
TOP_MARKERS=sys.argv[14]
EXPRESSION_FILTER=sys.argv[15]
MOTIF_FILTER=sys.argv[16]
TOP_REGULONS=sys.argv[17]

WARNING_MESSAGES = [] # will store warning messages to be displayed in the log


########## HELPER FUNCTIONS ##########

def path_is_dummy(path): # checks whether the user specified "null" as path to that file
        return path.startswith('.dummy_path')

def find_symbol_in_file(file_name, symbol):
    result = False
    try:
        output = subprocess.check_output(f"grep '{symbol}' {file_name}", shell=True)
        if output:
            result = True
    except subprocess.CalledProcessError as e:
        if e.returncode != 1:  # 1 means grep found no matches
            raise Exception(f"Command 'grep' (1) failed with the error: '{e}'!")
    return result


########## NUMERICAL PARAMETERS ##########

# CHECK: the parameters have correct values
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


########## GENERAL VERIFICATIONS ##########

# CHECK: if the GO file is null, then terms of interest is also null
if path_is_dummy(GO_FILE):
    if not path_is_dummy(TERMS_OF_INTEREST):
        raise Exception("When the GO file is set to 'null', the terms of interest file should also be set to null!")
    

# CHECK: if the TF info file is null, then the motif mapping step must be disabled
if path_is_dummy(INFO_TF):
    if DO_MOTIF_ANALYSIS == "true":
        raise Exception("When the TF info file is set to 'null', the motif mapping step should be disabled! (doMotifAnalysis = false)")


# CHECK: if the motif mapping file is null, then the motif mapping step must be disabled
if path_is_dummy(FEATURE_FILE_MOTIFS):
    if DO_MOTIF_ANALYSIS == "true":
        raise Exception("When the motif mapping file is set to 'null', the motif mapping step should be disabled! (doMotifAnalysis = false)")


# CHECK: equal number of files was provided for each file category, one per dataset
if len(EXPRESSION_MATRIX) != len(MARKERS_OUT):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of markers files ({len(MARKERS_OUT)})!")
if len(EXPRESSION_MATRIX) != len(CELLS2CLUSTERS):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of cell2clusters files ({len(CELLS2CLUSTERS)})!")
if len(EXPRESSION_MATRIX) != len(CLUSTER2IDENT):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of identities files ({len(CLUSTER2IDENT)})!")
if not path_is_dummy(GRNBOOST_OUT[0]) and len(EXPRESSION_MATRIX) != len(GRNBOOST_OUT):
    raise Exception(f"The number of expression matrix files ({len(EXPRESSION_MATRIX)}) is different from the number of GRNBoost output files ({len(GRNBOOST_OUT)})!")


# CHECK: input files, except the feature file motifs, are in plain text and not (g)zipped
try:
    all_dataset_file_names = EXPRESSION_MATRIX + MARKERS_OUT + CELLS2CLUSTERS + CLUSTER2IDENT
    all_file_names = all_dataset_file_names + GRNBOOST_OUT + [TF_LIST, TERMS_OF_INTEREST, INFO_TF, GO_FILE, ALIAS, ENRICHMENT_BACKGROUND]
    for file_name in all_file_names:
        if not path_is_dummy(file_name):
            output = subprocess.check_output(f'file "$(readlink -f {file_name})"', shell=True)
            if not "text" in output.strip().decode('utf-8'):
                raise Exception(f"The file '{file_name}' is compressed!")
        
    if not path_is_dummy(FEATURE_FILE_MOTIFS):
        output = subprocess.check_output(f'file "$(readlink -f {FEATURE_FILE_MOTIFS})"', shell=True)
        if not "gzip" in output.strip().decode('utf-8'):
            raise Exception(f"The file '{FEATURE_FILE_MOTIFS}' is not gzipped!")
except subprocess.CalledProcessError as e:
    raise Exception(f"Command 'file' (2) failed with the error '{e}'!")


########## DATASET-RELATED VERIFICATIONS ##########

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
    raise Exception("All input files must start with a dataset name followed by an underscore!")

# performs checks common to all input dataset files
def check_dataset_file(file_name: str):
    # CHECK: no quotes allowed in input files
    if find_symbol_in_file(file_name, '"'):
        raise Exception(f"Quotes (\") found in file'{file_name}'!")
        

    # CHECK: no trailing newlines present in the input files
    if find_symbol_in_file(file_name, '^$'):
        raise Exception(f"Empty lines found in '{file_name}'!")
    

    # CHECK: values in the first column of each data file must be unique
    # the only exception is the Seurat markers file, so it is not checked
    if file_name not in MARKERS_OUT:
        try:
            output = subprocess.check_output(f"cut -f1 {file_name} | sort -T . | uniq -d", shell=True)
            if output.strip():
                raise Exception(f"Duplicated indexes detected in file '{file_name}'!")
        except subprocess.CalledProcessError as e:
            raise Exception(f"Command 'cut' (3) failed with the error '{e}'!")
    
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
    if not path_is_dummy(GRNBOOST_OUT[0]) and not any(file_name.startswith(prefix) for file_name in GRNBOOST_OUT):
        raise Exception(f"No GRNBoost output file provided for the dataset {data_set}!")
    
    # retrieve the names of the dataset-related input files
    cells2clusters_file_name = [file_name for file_name in CELLS2CLUSTERS if file_name.startswith(prefix)][0]
    identities_file_name = [file_name for file_name in CLUSTER2IDENT if file_name.startswith(prefix)][0]
    matrix_file_name = [file_name for file_name in EXPRESSION_MATRIX if file_name.startswith(prefix)][0]
    markers_file_name = [file_name for file_name in MARKERS_OUT if file_name.startswith(prefix)][0]


    # CHECK: the right number of columns in the input files
    identifies_column_cnt = int(subprocess.check_output(f"head -n 1 {identities_file_name} | cut -f1- | wc -w", shell=True).strip())
    if identifies_column_cnt not in [2,3]:  # can have an optional third column with indices
        raise Exception(f"Incorrect number of columns in the file '{identities_file_name}': {identifies_column_cnt} instead of 2 (or 3)!")
    cells2clusters_column_cnt = int(subprocess.check_output(f"head -n 1 {cells2clusters_file_name} | cut -f1- | wc -w", shell=True).strip())
    if cells2clusters_column_cnt != 2:
        raise Exception(f"Incorrect number of columns in the file '{cells2clusters_file_name}': {cells2clusters_column_cnt} instead of 2!")
    markers_columns = subprocess.check_output(f"head -n 1 {markers_file_name}", shell=True).strip().decode('utf-8').split('\t')
    if not "p_val_adj" in markers_columns:
        raise Exception(f"The column 'p_val_adj' is missing in the file {markers_file_name}!")
    if not "cluster" in markers_columns:
        raise Exception(f"The column 'cluster' is missing in the file {markers_file_name}!")
    if not "gene" in markers_columns:
        raise Exception(f"The column 'gene' is missing in the file {markers_file_name}!")
        

    # CHECK: no underscores allowed in cluster names
    try:
        output = subprocess.check_output(f'cut -f2 {identities_file_name} | grep "_"', shell=True)
        if output:
            raise Exception(f"Underscores detected in cluster names of the '{data_set}' data set!")
    except subprocess.CalledProcessError as e:
        if e.returncode != 1:  # 1 means grep found no matches
            raise Exception(f"Command 'cut' (4) failed with the error: '{e}'!")
    

    # CHECK: cell ids correspond between the matrix and cells2clusters files
    try:
        # retrieve matrix cell ids from column names
        command = f"head -n1 {matrix_file_name}| tr '\t' '\n' | sort -T . -u"
        cells_from_matrix = set(subprocess.check_output(command, shell=True).decode().split())
    except subprocess.CalledProcessError as e:
        if e.returncode != 1:  # 1 means grep found no matches
            raise Exception(f"Command 'head' (5) failed with the error: '{e}'!")

    try:
        # retrieve cells2clusters cell ids from the first column
        command = f"cut -f1 {cells2clusters_file_name} | sort -T . -u"
        cells_from_cells_to_clusters = set(subprocess.check_output(command, shell=True).decode().split())
    except subprocess.CalledProcessError as e:
        if e.returncode != 1:  # 1 means grep found no matches
            raise Exception(f"Command 'cut' (6) failed with the error: '{e}'!")

    # check the intersection
    cells_in_common_count = len(cells_from_matrix & cells_from_cells_to_clusters)
    cells_in_matrix_count = len(cells_from_matrix)
    cells_in_cells_to_clusters_count = len(cells_from_cells_to_clusters)

    cell_id_example = ''  # random cell ID to show as example to the user
    if cells_in_matrix_count > 0:
        cell_id_example = next(iter(cells_from_matrix))
    elif cells_in_cells_to_clusters_count > 0:
        cell_id_example = next(iter(cells_from_cells_to_clusters))

    # all cell identities differ
    if cells_in_common_count == 0:
        raise Exception(f"The cell identifiers (e.g. {cell_id_example}) don't correspond between the matrix and cells2cluster file for the '{data_set}' data set!")
    # some cell identifiers differ
    elif not((cells_in_common_count == cells_in_matrix_count or
            cells_in_common_count == cells_in_matrix_count - 1) and   # in some cases the exported matrix is missing the first column (=index)
            cells_in_common_count == cells_in_cells_to_clusters_count):
        WARNING_MESSAGES.append(f"WARNING: The cell identifiers (e.g. {cell_id_example}) are not fully corresponding between the matrix and cells2cluster file for the '{data_set}' data set!")

    # CHECK: cluster ID and identities must correspond between cells2clusters and cluster2ident files
    try:
        output = subprocess.check_output(f'bash -c "comm -3 <(cut -f1 {identities_file_name} | sort -T . -u) <(cut -f2 {cells2clusters_file_name} | sort -T . -u)"', shell=True)
        if output:
            raise Exception(f"Cluster identities differ between cells2cluster and identities files for the '{data_set}' data set!")
    except subprocess.CalledProcessError as e:
        if e.returncode != 1:  # 1 means grep found no matches
            raise Exception(f"Command 'comm' (7) failed with the error: '{e}'!")


    # CHECK: only up-regulated markers should be present
    try:
        output = subprocess.check_output(f'cut -f3 {markers_file_name} | grep "-"', shell=True)
        if output:
            raise Exception(f"Down-regulated markers detected for the '{data_set}' data set!")
    except subprocess.CalledProcessError as e:
        if e.returncode != 1:  # 1 means grep found no matches
            raise Exception(f"Command 'cut' (8) failed with the error: '{e}'!")
    

    # CHECK: if the enrichment background is provided, then it must have at least one gene id in common with the expressed genes
    if not path_is_dummy(ENRICHMENT_BACKGROUND):
        try:
            output = subprocess.check_output(f'bash -c "comm -12 <(cut -f1 {matrix_file_name} | sort -T . -u) <(cut -f1 {ENRICHMENT_BACKGROUND} | sort -T . -u)"', shell=True)
            if not output:
                raise Exception(f"The enrichment background has no genes in common with the expressed genes in the '{data_set}' data set!")
        except subprocess.CalledProcessError as e:
            if e.returncode != 1:  # 1 means grep found no matches
                raise Exception(f"Command 'comm' (9) failed with the error: '{e}'!")


    # collect statistics
    stats_df.loc[data_set, 'cells'] = int(subprocess.check_output(f"head -n 1 {matrix_file_name} | cut -f1- | wc -w", shell=True).strip()) - 1
    stats_df.loc[data_set, 'genes'] = int(subprocess.check_output(f"cat {matrix_file_name} | wc -l", shell=True).strip()) - 1
    stats_df.loc[data_set, 'clusters'] = int(subprocess.check_output(f"cat {identities_file_name} | wc -l", shell=True).strip())
    stats_df.loc[data_set, 'tissues'] = int(subprocess.check_output(f"cut -f2 {identities_file_name} | sort -T . | uniq | wc -l", shell=True).strip())


# if the process is here, that means that all the tests passed, otherwise an exception is raized and the process is interrupted
print("")
print("== INPUT VALIDATION ==========================================")
print("Input files passed validation tests. Retrieved following data:")
print(stats_df)
for warning in WARNING_MESSAGES:  # 
    print(warning)
print("")
print("== INPUT FILES ===============================================")
print(f"Expression matrix file(s)  : {' / '.join(EXPRESSION_MATRIX)}")
print(f"Seurat markers file(s)     : {' / '.join(MARKERS_OUT)}")
print(f"Cells to clusters file(s)  : {' / '.join(CELLS2CLUSTERS)}")
print(f"Cluster identities file(s) : {' / '.join(CLUSTER2IDENT)}")
print(f"GRNBoost output file(s)    : {' / '.join(GRNBOOST_OUT) if not path_is_dummy(GRNBOOST_OUT[0]) else 'NOT PROVIDED'}")
print(f"Transcription factor file  : {TF_LIST}")
print(f"TF info file               : {INFO_TF}")
print(f"Gene aliases file          : {ALIAS}")
print(f"Motifs feature file        : {FEATURE_FILE_MOTIFS if not path_is_dummy(FEATURE_FILE_MOTIFS) else 'NOT PROVIDED'}")
print(f"GO file                    : {GO_FILE if not path_is_dummy(GO_FILE) else 'NOT PROVIDED'}")
print(f"Terms of interest file     : {TERMS_OF_INTEREST if not path_is_dummy(TERMS_OF_INTEREST) else 'NOT PROVIDED'}")
print(f"Enrichment background file : {ENRICHMENT_BACKGROUND if not path_is_dummy(ENRICHMENT_BACKGROUND) else 'NOT PROVIDED'}")
print("")
print("== MINI-EX PARAMETERS ========================================")
if not path_is_dummy(TERMS_OF_INTEREST):
    terms_of_interest_df = pd.read_csv(TERMS_OF_INTEREST, names=['term'])
    print(f"Terms of interest          : {' / '.join(terms_of_interest_df['term'])}")
print(f"doMotifAnalysis            : {DO_MOTIF_ANALYSIS}")
print(f"topMarkers                 : {TOP_MARKERS}")
print(f"expressionFilter           : {EXPRESSION_FILTER}")
print(f"motifFilter                : {MOTIF_FILTER}")
print(f"topRegulons                : {TOP_REGULONS}")
print("")