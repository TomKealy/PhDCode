function y = thresh(g,x)
% create a threshold graph, or test if g is a threshold graph
%
% thresh(g,x) -- create a threshold graph
% g is the graph to be created
% x is a vector of values in [0,1]
% we have an edge between i and j in g iff x(i) + x(j) >= 1
% 
% thresh(g) -- determine if g is a threshold graph; return true if it is.


% one argument case: check if the graph is a threshold graph
if nargin==1
    ds = sort(deg(g));
    y = dscheck(ds);
    return
end



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


function y = dscheck(ds)
% check a degree sequence to see if it comes from a threshold graph
% we assume ds is given to us sorted
if isempty(ds)
    y = true;
    return
end

if length(ds) == 1
    y = (ds == 0);
    return
end

if ds(1) == 0
    y = dscheck(ds(2:end));
    return
end

n = length(ds);
if (ds(n) == n-1)
    ds = ds(1:n-1)-1;
    y = dscheck(ds);
    return
end
y = false;
    


