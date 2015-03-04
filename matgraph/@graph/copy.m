function copy(g,h)
% copy(g,h) --- overwrite g with a copy of h
%
% This is different from writing g = h which makes g another view of the
% same graph h (so changes to h affect g). Rather, copy(g,h) makes g an
% identical, but independent copy of h.
%
% Note that g = graph(h) has a similar effect, but grabs a new slot in
% GRAPH_MAGIC. The syntax copy(g,h) is preferred when g is already a handle
% to a graph object.

global GRAPH_MAGIC

GRAPH_MAGIC.graphs{g.idx} = GRAPH_MAGIC.graphs{h.idx};

    