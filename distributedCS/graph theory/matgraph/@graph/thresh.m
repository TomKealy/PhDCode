function thresh(g,x)
% thresh(g,x) -- create a threshold graph
% g is the graph to be created
% x is a vector of values in [0,1]
% we have an edge between i and j in g iff x(i) + x(j) >= 1

% make sure x is a row-vector
x = x(:)';

n = length(x);

% Here is a MATLAB trick. We want a matrix whose ij entry is x(i)+x(j).
% To do this, we exponentiate x, multiple x'*x, and then logarithm. 

ex = exp(x);
M = log(ex'*ex);

A = M >= 1;
for k=1:n
	A(k,k) = 0;
end

set_matrix(g,A);