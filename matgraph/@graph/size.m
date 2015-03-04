function [sz,tz] = size(g)
% size(g) --- returns [nv,ne] for the graph

if nargout >= 2
    sz = nv(g);
    tz = ne(g);
else
    sz = [nv(g),ne(g)];
end