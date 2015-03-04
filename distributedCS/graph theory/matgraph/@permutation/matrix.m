function A = matrix(p)
% matrix(p) gives a permutation matrix version of p
n = length(p);
A = eye(n);
A = A(:,p.array);