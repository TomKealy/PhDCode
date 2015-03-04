function spy(g)
% spy(g) --- display the "spy" of g's adjacency matrix

global GRAPH_MAGIC
spy(GRAPH_MAGIC.graphs{g.idx}.array);