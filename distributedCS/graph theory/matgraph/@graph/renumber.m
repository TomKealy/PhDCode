function renumber(g,perm)
% renumber the vertices of a graph
% renumber(g,perm) --- renumber the vertices of a graph accoring to a
% permutation
% renumber(p,part) --- renumber vertices according to a partition.
%
% perm should be a permutation of 1 through n
% the graph's vertices are permutated so that the old vertex 1 is now
% vertex perm(1), etc.
% Note: perm may be either an array containing the permutation or a
% permutation object. 
%
% part should be a partition of [n]. We renumber so vertices in the same
% block of the partition are consecutive.


global GRAPH_MAGIC

n = nv(g);

if isa(perm,'partition')
    perm = partition2list(perm);
end

if ~isa(perm,'permutation')
    perm = permutation(perm);
end

if length(perm) ~= n
    error('Length of permutation does not match size of graph')
end

perm = inv(perm);
q = array(perm);

GRAPH_MAGIC.graphs{g.idx}.array = ...
    GRAPH_MAGIC.graphs{g.idx}.array(q,q);

make_logical(g);

if hasxy(g)
    xy = getxy(g);
    xy = xy(q,:);
    embed(g,xy);
end

if is_labeled(g)
    tmp = cell(n,1);
    for k=1:n
        tmp{k} = GRAPH_MAGIC.graphs{g.idx}.labels{perm(k)};
    end
    GRAPH_MAGIC.graphs{g.idx}.labels = tmp;
end


function list = partition2list(p)
pp = parts(p);
list = [];
for k=1:np(p)
    list = [list, pp{k}];
end
