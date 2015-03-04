function cycle(g,n)
% cycle(g,n) --- create a cycle graph on n vertices

path(g,n);
add(g,1,n);
embed(g);
