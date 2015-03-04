function free_all
% free_all --- delete all graphs from the GRAPH_MAGIC system

if ~graph_system_exists
    error('Graph system has not been intialized')
end

global GRAPH_MAGIC

for i=1:GRAPH_MAGIC.ngraphs
    GRAPH_MAGIC.graphs{i} = [];
    GRAPH_MAGIC.in_use(i) = 0;
end