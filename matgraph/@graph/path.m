function path(g,n)
% path(g,n) --- make g a path on n vertices


if hasxy(g)
    rmxy(g)
end

resize(g,0);
resize(g,n);

for k=1:n-1
    add(g,k,k+1);
end

x = 0:n-1;
y = zeros(1,n);
xy = [x',y'];

embed(g,xy);
