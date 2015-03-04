function m = laplacian(g,s)
% laplacian(g) --- get the Laplacian matrix of g
% equal to D-A where A is the adjacency matrix and D is a diagonal matrix
% of the degrees of the vertices.
% laplacian(g,'normalized') --- get the normalized Laplacian of g


n = nv(g);
d = deg(g);
m = -matrix(g);
for k=1:n
    m(k,k) = d(k);
end

if nargin>1
	S = zeros(n);
	for v=1:n
		if (d(v)>0)
			S(v,v) = 1/sqrt(d(v));
		else
			S(v,v) = 0;
		end
	end
	m = S*m*S;
end
