function full(g)
% full(g) --- convert internal storage for g to full

global GRAPH_MAGIC

GRAPH_MAGIC.graphs{g.idx}.array = ...
    full(GRAPH_MAGIC.graphs{g.idx}.array);