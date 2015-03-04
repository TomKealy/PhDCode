function out = get_label(g,v)
% get_label(g) or get_label(g,v) --- get vertex label(s)
% get_label(g) returns a cell array of the labels of the vertices
% get_label(g,v) returns v's label

global GRAPH_MAGIC;

if ~is_labeled(g)
    disp('graph is not labeled');
    out = {};
    if nargin==2
        out='';
    end
    return
end

if nargin==2
    out = GRAPH_MAGIC.graphs{g.idx}.labels{v};
    return
end

out = GRAPH_MAGIC.graphs{g.idx}.labels;