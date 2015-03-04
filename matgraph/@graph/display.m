function display(g)
% display(g) --- display basic information about a graph object

global GRAPH_MAGIC;

if (~GRAPH_MAGIC.in_use(g.idx))
    disp(['Invalid graph object (index ', int2str(g.idx),')']);
    return
end

n = nv(g);
m = ne(g);

if (isfull(g))
    stype = 'full';
else
    stype = 'sparse';
end


disp(['Graph with ', int2str(n), ' vertices and ', ...
    int2str(m),' edges (', stype, ')'])