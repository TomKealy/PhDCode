function circulant(g,n,k)
% circulant(g,n,k) --- overwrite g with an n,k circulant graph
% g has n vertices {1,2,...,n}. Each vertex has an edge to the next k in
% line. So this is a 2k regular graph.

resize(g,n);
clear_edges(g);

for u=1:n
    for uu=u+1:u+k
        v = mod(uu,n);
        if (v==0) 
            v = n;
        end
        add(g,u,v);
    end
end
