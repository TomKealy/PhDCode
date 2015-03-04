function nlist = neighbors(g,v)
% neighbors(g,v) --- neighbors of v as a list.
% Note: invoking g(v) has the same effect.

global GRAPH_MAGIC

n = nv(g);
if (v<0) || (v>n)
    nlist = [];
else
    nlist = find(GRAPH_MAGIC.graphs{g.idx}.array(v,:));
end
