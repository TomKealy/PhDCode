function nvec = subsref(g,S)
% subsref implements the notation g(v) and g(v,w)
% g(v) --- returns the neighbor set of v
% g(v,w) --- returns 1 if v and w are adjacent and 0 otherwise

global GRAPH_MAGIC

args = S.subs;
n = nv(g);

% single subscript case

if length(args)==1
    v = args(1);
    v = v{1};
    nvec = neighbors(g,v);
    return
end

% two subscript case

v = args(1);
v = v{1};
w = args(2);
w = w{1};

if (v<0) || (v>n) || (w<0) || (w>n) || (v==w)
    nvec = false;
else
    nvec = GRAPH_MAGIC.graphs{g.idx}.array(v,w);
    nvec = full(nvec);
end
