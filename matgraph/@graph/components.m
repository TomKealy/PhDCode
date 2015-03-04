function p = components(g)
% components(g) --- find the components of the graph g
% If g has n vertices, this returns a partition of the [n] based on the
% components of g.

n = nv(g);
indicator = zeros(n,1);
c = 0;

while (nnz(indicator)<n)
    c = c+1;
    % find first zero entry in indicator
    i = find(indicator==0);
    i = i(1);
    
    % get i's component
    ci = component(g,i);
    
    indicator(ci) = c;
end
%p = indicator;
p = ind2part(indicator);



function p = ind2part(ind)
% convert indicator vector to a partition

np = max(ind);
A = cell(np,1);
for k=1:np
    A{k} = find(ind == k);
end
p = partition(A);