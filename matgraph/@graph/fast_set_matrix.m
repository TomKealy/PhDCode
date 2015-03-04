function fast_set_matrix(g,A)
% fast_set_matrix(g,A) --- overwrite the adjacency matrix of g with A
% WARNING: No checks are done on the matrix!!! See set_matrix for a safer
% method. See also check_matrix.

global GRAPH_MAGIC
GRAPH_MAGIC.graphs{g.idx}.array = logical(A);