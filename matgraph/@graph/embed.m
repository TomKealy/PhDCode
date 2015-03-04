function embed(g,xy)
% embed --- create an embedding for a graph
% embed(g,xy) --- set the embedding to xy (an n-by-2 matrix)
% embed(g) --- default circulat embedding

global GRAPH_MAGIC;

n = nv(g);

if (nargin == 1)
    t = (0:n-1)*2*pi/n;
    x = n * cos(t)/6;
    y = n * sin(t)/6;
    xy = [x',y'];
end

[nr,nc] = size(xy);
if (nr ~= n) || (nc ~= 2)
    error('Embedding must be an n-by-2 matrix');
end

GRAPH_MAGIC.graphs{g.idx}.xy = xy;
