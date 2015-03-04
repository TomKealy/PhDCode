function p = swap(a,b)
% swap(a,b) -- return a permutation that swaps a and b, i.e., (a,b)

m = max(a,b);
row = 1:m;
row(a) = b;
row(b) = a;
p = permutation(row);
