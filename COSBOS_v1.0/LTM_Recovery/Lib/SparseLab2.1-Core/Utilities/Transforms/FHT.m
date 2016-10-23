function y = FHT(x)
% FHT: Computes the Hadamard transform of a signal
%  Usage:
%    y = FHT(x);
%  Inputs:
%    x      signal of dyadic length
%  Outputs:
%    y      Hadamard transform coefficients
%
% Description
%   FHT computes the Fast Hadamard transform of the signal x.
% See Also
%   IFHT

x = x(:);
n = length(x);

y = x;
t = zeros(size(x));  

odds = 1:2:n;
evens = 2:2:n;

for ii = 1:log2(n)
    t(odds) = y(odds) + y(evens);
    t(evens) = y(odds) - y(evens);
    y(1:(n/2)) = t(odds);
    y((n/2+1):n) = t(evens);
end

y = y ./ sqrt(n);
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
