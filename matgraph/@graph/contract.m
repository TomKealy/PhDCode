function contract(g,u,v)
% contract(g,u,v) --- contract v into u
% All neighbors of v are added to the neighborhood of u, 
% and then v is deleted. (When v is deleted, vertices that have higher
% index than v are renumbered.)

n = nv(g);

if (u<0) || (u>n) || (v<0) || (v>n) || (u==v) 
    return
end

nlist = neighbors(g,v);

for x = nlist;
    if (x ~= u)
        add(g,u,x)
    end
end

delete(g,v);
