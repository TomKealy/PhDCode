function bfstree(t,g,v)
% bfstree(t,g,v) --- create a breadth-first spanning tree of g
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

if (issparse(g))
    sparse(t)
end
resize(t,0);
resize(t,n);


q_init(n+1);

visited = zeros(n,1);
visited(v) = 1;
q_push(v);

while(q_size > 0)
    x = q_pop_front;
    nlist = neighbors(g,x);
    for w = nlist
        if (visited(w)==0)
            visited(w) = 1;
            add(t,x,w);
            q_push(w);
        end
    end
end


if hasxy(g)
    embed(t,getxy(g))
end

if is_labeled(g)
    copy_labels(t,g);
end

