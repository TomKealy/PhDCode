function clear_edges(g)
% clear_edges(g) --- delete all edges of g

global GRAPH_MAGIC

n = nv(g);

if (issparse(g))
    GRAPH_MAGIC.graphs{g.idx}.array = ...
        logical(sparse([],[],[],n,n,0));
else
    GRAPH_MAGIC.graphs{g.idx}.array = false(n);
end
