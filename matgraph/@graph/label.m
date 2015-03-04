function label(g,v,name)
% Assign labels to vertices of g
% With no arguments, we assign default names for vertices (string versions
% of the vertex numbers)
% label(g,v,name) gives vertex v the name in the string name.
% label(g,vnamelist) assigns names to vertices in the cell array vlist


global GRAPH_MAGIC;

n = nv(g);
i = g.idx;

if nargin == 1
    clear_labels(g);
    GRAPH_MAGIC.graphs{i}.labels = cell(n,1);
    for v=1:n
        GRAPH_MAGIC.graphs{i}.labels{v} = int2str(v);
    end
    return
end

if nargin == 2
    clear_labels(g);
    GRAPH_MAGIC.graphs{i}.labels = v;
    return    
end

if ~is_labeled(g)
    label(g);
end

if (v<0) || (v>n)
    error('No such vertex to label')
end

GRAPH_MAGIC.graphs{i}.labels{v} = name;