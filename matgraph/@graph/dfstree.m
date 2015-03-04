function dfstree(t,g,v)
% dfstree(t,g,v) --- create a depth-first spanning tree of g
% The tree is rooted a the vertex v (or vertex 1 if missing). If g is not
% connected, we generate a tree only for the component containing v;
% vertices in the other components are isolated vertices in t.



if nargin==2
    v = 1;
end

n = nv(g);
if (v<0) || (v>n)
    error('Seed vertex out of range')
end

copy(t,g);
clear_edges(t);

q_init(2*n+1);

visited = zeros(n,1);
tree_build(t,g,v,visited);
    
   
if is_labeled(g)
    copy_labels(t,g);
end


if hasxy(g)
    embed(t,getxy(g))
end

end


function visited = tree_build(t,g,v,visited)
% recursively build DFS tree

visited(v) = 1;

nlist = neighbors(g,v);

for w = nlist
    if (~visited(w))
        add(t,v,w);
        visited = tree_build(t,g,w,visited);
    end
end
        




end
