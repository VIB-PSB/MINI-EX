### adapted from aertslab/pySCENIC/src/pyscenic/cli/arboreto_with_multiprocessing.py

import sys
import time
import numpy
import pandas
from multiprocessing import Pool
import tqdm

from arboreto.algo import _prepare_input
from arboreto.core import SGBM_KWARGS, EARLY_STOP_WINDOW_LENGTH
from arboreto.core import target_gene_indices, infer_partial_network


tf_expression_matrix = sys.argv[1]
tf_list = sys.argv[2]
tg_expression_matrix = sys.argv[3]
tg_list = sys.argv[4]
num_workers = int(sys.argv[5])
output = sys.argv[6]
largest_int = int(sys.argv[7])
number_of_genes = int(sys.argv[8])
number_of_cells = int(sys.argv[9])


def run_infer_partial_network(target_gene_index):
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = ex_matrix[:, target_gene_index]

    n = infer_partial_network(
        regressor_type='GBM',
        regressor_kwargs=SGBM_KWARGS,
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=777)
    return( n )



if __name__ == '__main__':

    start_time = time.time()

    # read in TF matrix from binary numpy file
    tf_matrix = numpy.load(tf_expression_matrix)

    # read in the corresponding TF gene names
    with open(tf_list, 'r') as f:
        tf_names = [line.strip() for line in f]

    # check if the largest integer in the TG matrix is low enough to use int16 encoding, else use int32 
    encoding = numpy.uint16 if largest_int < 65536 else numpy.uint32

    # read in the expression matrix as a numpy array, assuming no header, transposing it
    ex_matrix = numpy.empty((number_of_cells, number_of_genes), dtype=encoding)
    with open(tg_expression_matrix, 'r') as f:
        for i in range(number_of_genes):
            line = f.readline().split('\t')
            ex_matrix[:, i] = numpy.array(line, dtype=encoding)

    # read in the corresponding gene names
    with open(tg_list, 'r') as f:
        gene_names = [line.strip() for line in f]
    
    end_time = time.time()
    print(f'Loaded expression matrix of {ex_matrix.shape[0]} cells and {ex_matrix.shape[1]} genes in {end_time - start_time} seconds...')
    print(f'Loaded {len(tf_names)} TFs...')

    ex_matrix, gene_names, _ = _prepare_input(ex_matrix, gene_names, None)

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
    