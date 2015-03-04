function [a,S] = alpha(g)
% [a,S] = alpha(g) --- indepencence number
% returns a, the independence number of the graph g, and an independent set
% S of maximum cardinality.
% Note: Uses linear programming and may work slowly (if at all) on large 
% graphs.
%
% REQUIRES THE OPTIMIZATION TOOLBOX

[n,m] = size(g);
M = incidence_matrix(g);
M = double(M');

c = -ones(1,n);
b = ones(m,1);

x = bintprog(c,M,b);

a = sum(x);
S = find(x)';
