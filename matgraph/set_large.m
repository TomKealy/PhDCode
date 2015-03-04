function sz = set_large(n)
% set_large(n) --- set the cut off size for large graphs. 
% set_large --- with no arguments, this returns the current threhold
% between "small" and "large".
%
% At construction, graphs on n or fewer vertices are considered "small" and
% are given full storage; otherwise, the graphs are "large" and given
% sparse storage.

if ~graph_system_exists
    error('Graph system not initialized')
end

global GRAPH_MAGIC

if nargin==0
    sz = GRAPH_MAGIC.large_size;
    return
end

if (n<1)
    n = 1000;
end

GRAPH_MAGIC.large_size = n;
sz = n;