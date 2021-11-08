### adapted from aertslab/pySCENIC/src/pyscenic/cli/arboreto_with_multiprocessing.py

import sys
import time
import pandas
from multiprocessing import Pool, cpu_count
import tqdm

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, _prepare_input
from arboreto.core import SGBM_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import to_tf_matrix, target_gene_indices, infer_partial_network



tfsList = sys.argv[1]
expressionMatrix = sys.argv[2]
num_workers = int(sys.argv[3])
output = sys.argv[4]



def run_infer_partial_network(target_gene_index):
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type='GBM',
        regressor_kwargs=SGBM_KWARGS,
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=777)
    return( n )



if __name__ == '__main__':

    start_time = time.time()
    ex_matrix = pandas.read_csv(expressionMatrix, sep='\t', header=0, index_col=0).T
    gene_names = ex_matrix.columns
    
    end_time = time.time()
    print(f'Loaded expression matrix of {ex_matrix.shape[0]} cells and {ex_matrix.shape[1]} genes in {end_time - start_time} seconds...')
    tf_names = load_tf_names(tfsList)
    print(f'Loaded {len(tf_names)} TFs...')

    ex_matrix, gene_names, tf_names = _prepare_input(ex_matrix, gene_names, tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(ex_matrix, gene_names, tf_names)
    start_time = time.time()

    with Pool(num_workers) as p:
        adjs = list(tqdm.tqdm(p.imap(run_infer_partial_network,
                                     target_gene_indices(gene_names, target_genes='all'),
                                     chunksize=1
                                     ),
                              total=len(gene_names)))

    adj = pandas.concat(adjs).sort_values(by='importance', ascending=False)

    end_time = time.time()
    print(f'Done in {end_time - start_time} seconds.')

    adj.to_csv(output, index=False, header=False, sep='\t')