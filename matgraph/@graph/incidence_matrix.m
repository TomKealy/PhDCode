function M = incidence_matrix(g,type)
% incidence_matrix(g) --- return the vertex/edge incidence matrix.
% The matrix returned is always sparse. 
%
% We return an nv-by-ne matrix whose ij entry is 1 if vertex i is an
% end point of edge j.
%
% Optional: incidence_matrix(g,'signed') which is the same matrix but one
% entry in each column is negative and one is negative.





[n,m] = size(g);
e = edges(g);

if nargin<2
	type = 'unsigned';  % default type
end

signed = false;

switch lower(type)
	case 'signed'
		signed = true;
	case 'unsigned'
		signed = false;
	otherwise
		error(['Unknown incidence matrix type:', type]);
end

i = [e(:,1);e(:,2)];
j = [1:m,1:m]';

if signed
	k = [ones(m,1);-ones(m,1)];
else
	k = ones(2*m,1);
end

M = sparse(i,j,k,n,m);
