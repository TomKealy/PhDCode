function complement(g,h)
% graph complement
% complement(g) --- overwrite a graph with its own complement
% complement(g,h) --- overwrite g with the complement of h

if nargin>1
    copy(g,h)
end

% complement only makes sense for full storage
if (issparse(g))
    full(g);
end

n = nv(g);
i = g.idx;

fast_set_matrix(g, ones(n) - eye(n) - matrix(g))