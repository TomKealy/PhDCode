function R = resistance(g)
% R = resistance(g) --- calculate resistances between vertices
% g is a connected, undirected weighted graph. Returns a matrix R whose 
% (i,j)-entry is the effective resistance between vertices i and j 
% where we assume each edge represents a resistance equal to the weight 
% of the edge.
%
% Author: Somit Gupta

if(isconnected(g)) %to test for connected graphs as the formulae works only for connected graphs
	l= laplacian(g); % gets the laplacian of the graph
	m=nv(g); %no of vertices
	% R=[0];
	lp=pinv(l);

	d=diag(lp);
	rowMatrix=d*ones(1,m);
	R = rowMatrix+rowMatrix'-2*lp;
else
    error('The graph is not connected')
end
