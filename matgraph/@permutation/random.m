function q = random(p)
% random(p) --- return a random reordering of the elements of p
a = randperm(length(p));
q = permutation(a);