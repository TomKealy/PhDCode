function nt = nsptrees(g)
% nsptrees(g) --- number of spanning trees of g.
% This computes the number of spanning trees of g using the matrix-tree
% theorem.

A = diag(deg(g)) - matrix(g);
n = nv(g)-1;
nt = det(A(1:n,1:n));