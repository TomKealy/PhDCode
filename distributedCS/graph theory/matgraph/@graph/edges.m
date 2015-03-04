function elist = edges(g)
% edges(g) --- list the edges in g as a 2-column matrix

global GRAPH_MAGIC
[i,j] = find(triu(GRAPH_MAGIC.graphs{g.idx}.array));

elist = [i,j];