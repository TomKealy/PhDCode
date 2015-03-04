function induce(g,arg1,arg2)
% induce --- create induced subgraphs
% induce(g,vset) --- g overwritten with induced subgraph on vset
% induce(g,h,vset) --- g overwritten with induced subgraph of h on vertex
%                      set vset


global GRAPH_MAGIC

if nargin == 2
    vset = arg1;
    if hasxy(g)
        GRAPH_MAGIC.graphs{g.idx}.xy = ...
            GRAPH_MAGIC.graphs{g.idx}.xy(vset,:);
    end
    if is_labeled(g)
        GRAPH_MAGIC.graphs{g.idx}.labels =  ...
            GRAPH_MAGIC.graphs{g.idx}.labels(vset);
    end
    GRAPH_MAGIC.graphs{g.idx}.array = ...
        GRAPH_MAGIC.graphs{g.idx}.array(vset,vset);
    return
end

h = arg1;
vset = arg2;

GRAPH_MAGIC.graphs{g.idx}.array = ...
    GRAPH_MAGIC.graphs{h.idx}.array(vset,vset);

if (hasxy(h))
    xy = getxy(h);
    xy = xy(vset,:);
    embed(g,xy);
end

if is_labeled(g)
    label(h);
    GRAPH_MAGIC.graphs{g.idx}.labels = ...
        GRAPH_MAGIC.graphs{h.idx}.labels(vset);
end