function cells = parts(p)
% parts(p) --- returns the parts of p in a cell array

[m,n] = size(p.array);

cells = cell(1,m);

for k=1:m
    cells{k} = find(p.array(m+1-k,:));
end
