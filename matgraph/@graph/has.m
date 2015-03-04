function yn = has(g,u,v)
% has --- check if the graph has a particular vertex or edge
% has(g,u) --- check if g contains the vertex u
% has(g,u,v) --- check if g contains the edge uv

global GRAPH_MAGIC;
n = nv(g);

if (nargin==2)
    yn = (u>0) & (u<=n);
else
    if (u<0) || (u>n) || (v<0) || (v>n) || (u==v)
        yn = 0;
    else
        yn = GRAPH_MAGIC.graphs{g.idx}.array(u,v);
    end
end
