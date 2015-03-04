function max_graphs = max_available
% num_available --- maximum number of graphs we can hold in this system

if ~graph_system_exists
    error('Graph system has not been initialized.')
end

global GRAPH_MAGIC;
max_graphs = GRAPH_MAGIC.ngraphs;