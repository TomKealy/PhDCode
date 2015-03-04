function elist = across(g,X,Y)
% across(g,X,Y) --- find all edges with one end in X and other end in Y.
% We require X and Y to be disjoint. If Y is omitted, then Y = V(g)-X.
% We do not require X U Y = V(g).

A = matrix(g);
X = X(:);
n = nv(g);

if nargin==2
    Y = setdiff(1:n,X);
end

Y = Y(:);

B = A(X,Y);
[i,j] = find(B);

ii = X(i);
jj = Y(j);
% make sure these are column vectors
ii = ii(:);
jj = jj(:);

elist = [ii,jj];

