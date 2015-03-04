function yn = isconnected(g)
% isconnected(g) --- test if g is a connected graph

global GRAPH_MAGIC

n = nv(g);
if (n==0) 
    yn = 1;
    return
end

x = zeros(n,1);
x(1) = 1;
while(1)
    y = x;
    x = GRAPH_MAGIC.graphs{g.idx}.array*x + x;
    x = double(x>0);
    if (x==y)
        break
    end
end

yn = sum(x)==n;
