function match_join(g,h1,h2)
% match_join(g,h1,h2) --- overwrite g with a graph formed from disjoint
% graphs g1,g2 with a matching added between them. The graphs h1 and h2
% must have the same number of vertices. 

if (nv(h1) ~= nv(h2))
    error('In match_join(g,h1,h2) graphs h1 and h2 must have the same number of vertices');
end

A = matrix(h1);
B = matrix(h2);
n = nv(h1);

set_matrix(g, [A, eye(n); eye(n), B]);