function cycle(g,n)
% cycle(g,n) --- create a cycle graph on n vertices

global GRAPH_MAGIC

path(g,n);
add(g,1,n);
embed(g);