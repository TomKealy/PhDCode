function m = ne(g,h)
% ne(g) --- number of edges in g or ne(g,h) --- check for inequality

global GRAPH_MAGIC

if (nargin == 2)
    if any(size(g) ~= size(h))
        m = 1;
        return
    end
    D = GRAPH_MAGIC.graphs{g.idx}.array - GRAPH_MAGIC.graphs{h.idx}.array;
    m = nnz(D)>0;
    return
end

m = nnz(GRAPH_MAGIC.graphs{g.idx}.array)/2;