function p = permutation(n)
% constructor for the permutation class
% p = permutation(n) --- create a new (identity) permutation on 1:n
% p = permutation(vec) --- create a new permutation based on the array vec
% (which should be a permutation of 1:n).

if nargin==0
    n = 1;
end

bad = 'Input invalid: must be a permutation of 1:n or a single integer';
[r,c] = size(n);

if (r>1) && (c>1)
    error(bad);
end

if (r==1) && (c==1)
    x = 1:n;
else
    x = n;
end


x = x(:)';
n = length(x);


y = sort(x);
if any(y ~= 1:n)
    error(bad);
end

p.array = x;
p = class(p,'permutation');
