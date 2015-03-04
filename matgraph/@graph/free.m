function free(g)
% free(g) --- free the graph from the system

if ~graph_system_exists
    error('Graph system has not been initialized')
end

global GRAPH_MAGIC
i = g.idx;

GRAPH_MAGIC.in_use(i) = 0;
GRAPH_MAGIC.graphs{i} = [];