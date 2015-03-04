function clear_labels(g)
% clear_labels(g) --- delete all labels in g

global GRAPH_MAGIC;
% GRAPH_MAGIC.graphs{g.idx}.labels = cell(0);

if isfield(GRAPH_MAGIC.graphs{g.idx},'labels')
	GRAPH_MAGIC.graphs{g.idx} = ...
		rmfield(GRAPH_MAGIC.graphs{g.idx},'labels');
end
