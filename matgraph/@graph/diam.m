function d = diam(g)
% diam(g) --- return the diameter of g: the length of a longest shortest
% path in g, or Inf if g is not connected.

d = max(max(dist(g)));