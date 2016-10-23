function alpha = GenRandVec(n,p)
% GenRandDict: Generates a random vector in R^n,
% with p-norm equal to C = (log(n))^(1/p).

C = log(n).^(1/p);
t = C .* ((1:n).*log(n)) .^ (-1/p);
tt = t(randperm(n));
alpha = (sign(randn([1 n])) .* tt).';

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
