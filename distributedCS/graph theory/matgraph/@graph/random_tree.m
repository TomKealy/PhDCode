function random_tree(t,n)
% random_tree(t,n) --- overwrite t with a random tree on n vertices

resize(t,n);
clear_edges(t);

for v=2:n
    u = floor(rand(1)*(v-1))+1;
    add(t,u,v);
end

p = randperm(n);
renumber(t,p);