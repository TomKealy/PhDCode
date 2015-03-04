function rmxy(g)
% rmxy(g) --- erase g's embedding

global GRAPH_MAGIC

if hasxy(g)
    GRAPH_MAGIC.graphs{g.idx} = rmfield(GRAPH_MAGIC.graphs{g.idx},'xy');
end