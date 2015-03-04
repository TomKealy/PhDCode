function set_matrix(g,A)
% set_matrix(g,A) --- set g to be the graph specificed in the matrix A.
% The matrix is scanned by check_matrix to be sure it is a valid adjacency
% matrix. If it is not, an error is triggered and the graph is left
% unchanged. 

if (~check_matrix(A))
    error('Input matrix is not a valid adjacency matrix')
end

fast_set_matrix(g,logical(A));

