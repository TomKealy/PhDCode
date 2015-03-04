function A = matrix(g)
% matrix(g) --- return (a copy of) the adjacency matrix of g
global GRAPH_MAGIC
A = GRAPH_MAGIC.graphs{g.idx}.array;