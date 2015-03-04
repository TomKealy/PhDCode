function r = mtimes(p,q)
% p*q is the composition of p and q
np = length(p);
nq = length(q);

if np < nq
    p.array = [p.array, np+1:nq];
    np = nq;
end

if np > nq
    q.array = [q.array, nq+1:np];
    % nq = np;
end

n = np;

a = zeros(1,n);

for k=1:n
    a(k) = p.array(q.array(k));
end

r = permutation(a);
