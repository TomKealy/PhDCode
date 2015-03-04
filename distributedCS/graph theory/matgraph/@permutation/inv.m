function q = inv(p)
% inv(p) returns the inverse permutation of p
a = array(p);
[y,b] = sort(a);
q = permutation(b);