function n = nv(g) 
% nv(g) --- number of vertices in g

global GRAPH_MAGIC
n = length(GRAPH_MAGIC.graphs{g.idx}.array);
