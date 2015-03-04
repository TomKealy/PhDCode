function yn = is_labeled(g)
% is_labeled(g) --- determine if there are labels on vertices.
% By definition, an empty graph is NOT labeled.

global GRAPH_MAGIC

if (isfield(GRAPH_MAGIC.graphs{g.idx},'labels'))
    yn = length(GRAPH_MAGIC.graphs{g.idx}.labels) > 0;
else
    yn = false;
end