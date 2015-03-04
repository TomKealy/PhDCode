function yn = hasxy(g)
% hasxy(g) --- determine if an embedding has been created for g

global GRAPH_MAGIC;

n = nv(g);

if isfield(GRAPH_MAGIC.graphs{g.idx},'xy')
    % the field exsts. now is is the right size?
    [nr,nc] = size(GRAPH_MAGIC.graphs{g.idx}.xy);
    if (nr == n) && (nc == 2)
        yn = 1;
        return
    else
        yn = 0;
        GRAPH_MAGIC.graphs{g.idx} = ...
            rmfield(GRAPH_MAGIC.graphs{g.idx},'xy');
        return
    end
end
yn = 0;
