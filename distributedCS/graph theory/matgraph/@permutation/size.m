function s = size(p)
% size(p) returns the number of elements permuted and the number of cycles.
s1 = length(p.array);
c = cycles(p);
s2 = length(c);
s = [s1,s2];