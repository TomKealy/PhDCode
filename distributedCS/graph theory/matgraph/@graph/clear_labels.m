function clear_labels(g)
% clearl_labels(g) --- delete all labels in g

global GRAPH_MAGIC;
GRAPH_MAGIC.graphs{g.idx}.labels = cell(0);
rmfield(GRAPH_MAGIC.graphs{g.idx},'labels');