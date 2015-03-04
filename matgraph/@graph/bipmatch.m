function elist = bipmatch(g,X,Y) 
% bipmatch --- maximum matching in a bipartite graph
% This may be called either as bipmatch(g,X,Y) or bipmatch(g,p).
% In either case, g is any graph.
% X,Y are disjoint sets of vertices of g or p is a partition of a subset of
% V(g) into two parts.
% In either case, this returns a maximum matching on the bipartite 
% (sub)graph of G with bipartition (X,Y).


if nargin==2
    c = parts(X);
    if length(c) ~= 2
        error('partition must have two parts');
    end
    X = c{1};
    Y = c{2};
else
    X = unique(X);
    Y = unique(Y);
end

if length(intersect(X,Y)) > 0
    error('The sets X and Y must be disjoint')
end

n = nv(g);

Z = union(X,Y);

if any(Z<1) || any(Z>n)
    error('Some vertices out of range')
end

A = matrix(g);
A = A(X,Y);

p = dmperm(A); % one of Matlab's best kept secrets :-)

yidx = find(p);
xidx = p(p>0);

x = X(xidx);
y = Y(yidx);

elist = [ x(:), y(:)];
elist = sortrows(elist);
