function fgraph(g,n,func)
% fgraph(g,n,func) -- create a graph based on a function.
% n is the number of vertices 
% func is the name of a function that is called with three arguments:
% The numbrer of vertices n, and two vertices u,v: func(n,u,v)
%
% We create a graph g that has n vertices and an edge from u to v iff
% func(n,u,v) is nonzero. The function is supplied by the user and it is the
% user's responsibility that this function be symmetric in u and v.

resize(g,0)
resize(g,n)

A = zeros(n,n);

for u=1:n-1
    for v=(u+1):n
        A(u,v) = feval(func,n,u,v);
    end
end
A = A+A';
A = A~=0;
set_matrix(g,A);