function s = size(p)
% size(p) --- returns [nv,np]: number of vertices, parts.
[m,n] = size(p.array);
s = [n,m];
