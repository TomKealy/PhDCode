function randxy(g)
% randxy(g) --- give g a random embedding in the plane

n = nv(g);
xy = sqrt(n)*randn(n,2)/2;
embed(g,xy);