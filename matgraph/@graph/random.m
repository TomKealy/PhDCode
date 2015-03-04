function random(g,arg1,arg2)
% random --- overwrite g with a random graph
% random(g,n,p) --- overwrite g with an Erdos-Renyi G(n,p) graph
% random(g,p) --- overwrite g with an Erdos-Reniy G(n,p) graph where n is
%                 the number of vertices already in g.
% random(g) --- equivalent to random(g,1/2)



% case: random(g,n,p)

if nargin==3
    n = arg1;
    p = arg2;
    
    A = triu(rand(n,n)<=p);
    A = A+A';
    for k=1:n
        A(k,k)=0;
    end
    fast_set_matrix(g,A);
    return
end

% case: random(g,p)

if nargin==2
    n = nv(g);
    p = arg1;
    random(g,n,p);
    return
end

% last case: random(g)

random(g,0.5);
