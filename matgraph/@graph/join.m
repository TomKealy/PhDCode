function join(g,h1,h2)
% join(g,h1,h2) --- overwrite g with the join of h1 and h2

global GRAPH_MAGIC

n1 = nv(h1);
n2 = nv(h2);

z = true(n1,n2);

rmxy(g);


GRAPH_MAGIC.graphs{g.idx}.array = ...
    [GRAPH_MAGIC.graphs{h1.idx}.array, z; ...
     z', GRAPH_MAGIC.graphs{h2.idx}.array
    ];


if (hasxy(h1) && hasxy(h2))
    xy = [getxy(h1);getxy(h2)];
    embed(g,xy);
else
	rmxy(g)
end
