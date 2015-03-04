function intersect(g,h1,h2)
% intersect(g,h1,h2) --- g is set to the intersection of h1 and h2.
% The vertex set of g is intersection of the vertex sets of h1 and h2, and
% the edge set of g is the intersection of the edge sets of h1 and h2.

global GRAPH_MAGIC


n1 = nv(h1);
n2 = nv(h2);

g1 = graph;
g2 = graph;

copy(g1,h1);
copy(g2,h2);

if n1 > n2
    resize(g1,n2)
else
    resize(g2,n1)
end

GRAPH_MAGIC.graphs{g.idx}.array = ...
    (matrix(g1) & matrix(g2));

free(g1);
free(g2);
