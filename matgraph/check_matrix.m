function ok = check_matrix(A)
% check_matrix(A) --- check if A is a valid adjacency matrix for a graph

ok = 0;

% make sure matrix is equare
[nr,nc] = size(A);
if (nr ~= nc) 
    return
end

% see if the matrix is zero-one
B = (A>0); % convert to logical (so 0,1).
if (nnz(A-B)>0)
    return
end

% check if symmetric
if (nnz(A-A')>0)
    return
end

% check the diagonal
if (any(diag(A)>0)) 
    return
end

% all tests passed
ok = 1;