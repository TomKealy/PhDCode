function sparse(g)
% sparse(g) --- convert internal storage for g to sparse

global GRAPH_MAGIC

GRAPH_MAGIC.graphs{g.idx}.array = ...
    sparse(GRAPH_MAGIC.graphs{g.idx}.array);