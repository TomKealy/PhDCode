function yn = issparse(g)
% issparse(g) --- check if g's adjacency matrix is sparse
global GRAPH_MAGIC
yn = issparse(GRAPH_MAGIC.graphs{g.idx}.array);