"""
Reads the TF expression matrix into memory as a numpy array, and saves it as a binary numpy file
for quick access.
"""

import sys,numpy

matrix_file = sys.argv[1]
output_file = sys.argv[2]
largest_int = int(sys.argv[3])
number_of_genes = int(sys.argv[4])
number_of_cells = int(sys.argv[5])

# check if the largest integer is low enough to use int16 encoding, else use int32 
encoding = numpy.uint16 if largest_int < 65536 else numpy.uint32

# read in the expression matrix as a numpy array, assuming no header, transposing it (since GRNBoost2 expects cells as rows and genes as columns)
expression_matrix = numpy.empty((number_of_cells, number_of_genes), dtype=encoding)
with open(matrix_file, 'r') as f:
    for i in range(number_of_genes):
        line = f.readline().split('\t')
        expression_matrix[:, i] = numpy.array(line, dtype=encoding)

# save the transposed matrix as binary numpy array
numpy.save(output_file, expression_matrix)
