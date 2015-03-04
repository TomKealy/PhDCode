function p = bipartition(g)
% part = bipartition(g) --- find bipartition of a bipartite graph
% If the graph is bipartite, this function finds a bipartition X,Y of the
% vertex set; these are returned in a partition object.
% If the graph is not bipartite, an empty partition is returned.

n = nv(g);

if n < 2
    p = partition(n);
    return
end

if ne(g) == 0
    p1 = 1:2:n;
    p2 = 2:2:n;
    p = partition({p1,p2});
    return
end

c = components(g);
x = zeros(n,1);
y = x;

nc = np(c);
c = parts(c);  % convert to cell array

% set up the X-seed vector
for k=1:nc
    j = min(c{k});
    x(j) = 1;
end

A = matrix(g);

last_sum = 0;
while sum(x) > last_sum
    last_sum  = sum(x) ;
    xx = A*(A*x)+x;
    x = double(xx>0);
end


yy = A*x;
y = double(yy>0);

% kludge to deal with isolated vertices
%idx = find((x==0)&(y==0));
%if ~isempty(idx)
%    x(idx) = 1;
%end




if sum(x) + sum(y) ~= n
    p = partition(0);
else
    part = cell(1,2);
    part{1} = find(x)';
    part{2} = find(y)';
    p = partition(part);
end