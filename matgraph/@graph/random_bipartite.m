function random_bipartite(g,n,m,p)
% random_bipartite(g,n,m,p) --- create a random bipartite graph
% the graph g is overwritten with a bipartite grapn with part sizes n and m
% with edge probability p. If p is missing, we take p = 1/2.

if nargin == 3
    p = 0.5;
end


A = rand(n,m) <= p;
fast_set_matrix(g,  [ zeros(n,n), A; A', zeros(m,m) ] );


rmxy(g);

xy1 = [ (1:n)-n/2; -ones(1,n)*log(n+m+1) ]';
xy2 = [ (1:m)-m/2; ones(1,m)*log(n+m+1) ]';
xy = [xy1;xy2];
embed(g,xy);


clear_labels(g);
