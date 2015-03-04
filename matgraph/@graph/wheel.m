function wheel(g,n)
% wheel(g,n) --- overwrite g with a wheel graph on n vertices
% if n is omitted, we use the number of vertices in g

if (nargin == 1)
    n = nv(g);
end

cycle(g,n-1);
resize(g,n);
for k=1:n-1
    add(g,k,n);
end