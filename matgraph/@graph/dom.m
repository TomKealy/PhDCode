function [d,S] = dom(g)
% [d,S] = dom(g) --- domination number of g
% The domination number is the smallest size of a set of vertices S with
% the property that every vertex is either in S or else distance 1 from
% some vertex in S. 
% Returns d, the domination number, and S a smallest dominating set.
% WARNING: It is NP hard to find the domination number of a graph, so this
% will work slowly (if at all) for large graphs. 
%
% REQUIRES THE OPTIMIZATION TOOLBOX

n = nv(g);
A = matrix(g) + eye(n);
f = ones(1,n);
b = ones(n,1);
[x,d] = bintprog(f,-A,-b);
S = find(x);
