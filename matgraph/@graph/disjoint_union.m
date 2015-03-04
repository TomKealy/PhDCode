function disjoint_union(g,h,k)
% disjoint_union(g,h,k) --- overwrite g with the disjoint union of the
% graphs h and k. 
% If h and k both have embeddings, those embeddings are combined to form an
% embedding for g.

global GRAPH_MAGIC

n1 = nv(h);
n2 = nv(k);

if (issparse(h) || issparse(k))
    z = logical(sparse([],[],[],n1,n2,0));
else
    z = false(n1,n2);
end

rmxy(g);

GRAPH_MAGIC.graphs{g.idx}.array = ...
    [GRAPH_MAGIC.graphs{h.idx}.array, z; ...
     z', GRAPH_MAGIC.graphs{k.idx}.array
    ];
make_logical(g);

if (hasxy(h) && hasxy(k))
    xy = [getxy(h);getxy(k)];
    embed(g,xy);
end
