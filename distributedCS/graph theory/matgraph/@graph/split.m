function part = split(g)
% part = split(g) --- partition the vertices of g into two subsets
% The graph is partitioned based on the eigenvector associated with the 2nd
% smallest eigenvalue of the Laplacian of g. 

L = laplacian(g);

[V,d] = eig(L);
v2 = V(:,2);

% if (sum(v2)<0) 
%     v2 = -v2;
% end

part = cell(1,2);

part{1} = find(v2 > 0);
part{2} = find(v2 <= 0);
part = partition(part);
