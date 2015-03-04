function v = eig(g)
% eig(g) -- compute the eigenvalues of this graph
% This returns the eigenvalues of the adjacency matrix of this graph.
% It's equivalent to eig(double(matrix(g)))

v = eig(double(matrix(g)));