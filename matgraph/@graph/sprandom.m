function sprandom(g,n,p)
% sprandom(g,n,p) --- create a sparse random graph
% This overwrites g with a sparse random graph on n vertices with edge
% probability p. 

global GRAPH_MAGIC

A = triu(sprand(n,n,p));
A = A + A';
for k=1:n
    A(k,k)=0;
end
A = logical(A>0);
GRAPH_MAGIC.graphs{g.idx}.array = A;
rmxy(g);