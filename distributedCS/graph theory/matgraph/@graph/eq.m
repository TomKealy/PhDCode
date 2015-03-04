function yn = eq(g,h)
% yn = eq(g,h) --- check if two graphs are equal (identical)

global GRAPH_MAGIC

if (nargin == 2)
    if any(size(g) ~= size(h))
        yn = 0;
        return
    end
    D = GRAPH_MAGIC.graphs{g.idx}.array - GRAPH_MAGIC.graphs{h.idx}.array;
    yn = nnz(D)==0;
    return
end