function paley(g,n)
% paley(g,n) --- create a Paley graph with n vertices
% In this graph, there is an edge from i to j iff i-j is a square in Z_n.
% Hence, if -1 is not a square in Z_n, an error is generated.

if (n <= 1) 
    error('n must be at least 2');
end

squares = unique(mod((1:n-1).^2,n));

% check if -1 is a square
if ~ismember(n-1,squares)
    error('-1 is not a square in Z_n');
end
resize(g,0);
resize(g,n);
rmxy(g);
for v = 1:n
    nlist = mod(v+squares,n);
    idx = find(nlist == 0);
    if length(idx)>0
        nlist(idx) = n;
    end
    for w = nlist
        add(g,v,w)
    end
end
