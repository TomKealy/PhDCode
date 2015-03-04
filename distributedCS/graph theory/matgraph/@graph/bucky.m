function bucky(g)
% bucky(g) --- overwrite g with the Buckyball graph (and give a nice
% embedding).

[A,xyz] = bucky;
set_matrix(g,full(A));
embed(g,2*xyz(:,1:2));