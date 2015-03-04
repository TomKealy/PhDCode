function p = chromatic_poly(g)
% chromatic_poly(g) --- find the chromatic polynomial of g
% Warning: This algorithm is extremely slow and unusable except for very
% small graphs.

n = nv(g);
m = ne(g);

if m==0
    p = zeros(1,n+1);
    p(1) = 1;
    return
end

elist = edges(g);
u = elist(1,1);
v = elist(1,2);

h = graph;
copy(h,g);

delete(h,u,v);
p1 = chromatic_poly(h);
contract(h,u,v);
p2 = chromatic_poly(h);

p = p1 - [0,p2];

free(h);