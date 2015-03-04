function nfree = num_available
% num_available --- number of free slots in the graph system

if ~graph_system_exists
    error('Graph system has not been initialized.')
end

global GRAPH_MAGIC;
nfree = GRAPH_MAGIC.ngraphs - sum(GRAPH_MAGIC.in_use);