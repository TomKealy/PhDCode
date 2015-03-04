function m = laplacian(g)
% laplacian(g) --- get the Laplacian matrix of g
% equal to D-A where A is the adjacency matrix and D is a diagonal matrix
% of the degrees of the vertices.


n = nv(g);
d = deg(g);
m = -matrix(g);
for k=1:n
    m(k,k) = d(k);
end