function DD = mdsxy(g,D)
% mdsxy(g) -- create an embedding based on multidimensional scaling
% 
% This may also be invoked as mdsxy(g,D) where D is a distance matrix
% (i.e., the i,j entry of D is the desired distance between i and j). 
% The default is the graph-theoretic distance from i to j.

% if the user does not 
if nargin == 1
	D = dist(g);
end

% any distances > cutoff  are truncated to cutoff.

cutoff = 10;
[i,j] = find(D>cutoff);
m = length(i);
for k=1:m
	D( i(k), j(k) ) = cutoff;
end



% rescale the D matrix so that the smallest positive entry is 1

vals = D(:);
vals = sort(vals);
vals = vals(vals>0);
minval = vals(1);
D = D/minval;


% find embedding vectors
D = D.^2;
n = nv(g);

d1 = sum(D)/n;
d2 = sum(d1)/n;

u = ones(1,n);

D1 = u'*d1;
D2 = d2 * ones(n,n);

Dstar = (-1/2) * (D - D1 - D1' + D2);

[v,d] = eig(Dstar);
d = real(d);
v = real(v);
dvec = (diag(d));
[dvec,idx] = sort(-dvec);

a = idx(1);
b = idx(2);

Lambda = d([a,b],[a,b]);
sqrt(Lambda);
xy = v(:,[a,b])*sqrt(Lambda);


embed(g,xy)

if nargout > 0
	DD = Dstar;
end
