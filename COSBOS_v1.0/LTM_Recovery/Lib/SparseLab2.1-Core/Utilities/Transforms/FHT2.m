function y = FHT2(x)
% FHT2: Computes the 2-D Hadamard transform of a signal
%  Usage:
%    y = FHT2(x);
%  Inputs:
%    x      image of dyadic dimensions
%  Outputs:
%    y      Hadamard transform coefficients
%
% Description
%   FHT2 computes the Fast 2-D Hadamard transform of the image x.
% See Also
%   FHT

[m,n] = size(x);
y = zeros(m,n);

for jj = 1:n
    y(:,jj) = FHT(x(:,jj));
end

for jj = 1:m
    y(jj,:) = FHT(y(jj,:))';
end
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
