function ds = vertex_degree_sequence(g,v)
% find the degrees of the neighbors of a vertex
vv = neighbors(g,v);
ds = [];
for w=vv
	ds = [ds,deg(g,w)];
end
ds = -sort(-ds);
