function yn = isempty(p)
% isempty(p) --- check if the permutation p is of size 0
s = size(p);
yn = s(1)==0;
