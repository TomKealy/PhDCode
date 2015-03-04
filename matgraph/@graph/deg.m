function d = deg(g,v) 
% deg: degree of a vertex or degree sequence
% deg(g) --- degree sequence of g
% deg(g,v) --- degree of v in g

global GRAPH_MAGIC


if (nargin==1)
    d = full(sum(GRAPH_MAGIC.graphs{g.idx}.array));
else
    d = full(sum(GRAPH_MAGIC.graphs{g.idx}.array(v,:)));
end