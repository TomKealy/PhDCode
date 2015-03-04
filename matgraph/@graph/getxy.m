function z = getxy(g)
% getxy(g) --- give g's embedding (or [] if g doesn't have one)


global GRAPH_MAGIC

if (hasxy(g))
    z = GRAPH_MAGIC.graphs{g.idx}.xy;
else
    z = [];
end