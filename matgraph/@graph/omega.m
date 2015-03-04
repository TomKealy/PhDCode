function [w,S] = omega(g)
% [w,S] = omega(g) --- clique number
% Returns w, the size of a largest clique, and S, a largest clique.
% Note: Uses integer programming so may run slowly (if at all) on larger
% graphs. 
%
% REQUIRES THE OPTIMIZATION TOOLBOX

% Note: Our implementation is simply alpha(-g).

h = graph;
copy(h,g);
complement(h);
[w,S] = alpha(h);
free(h);
