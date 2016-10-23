function y = RST(x)

y = fftshift(fft(fftshift(x)))./ sqrt(prod(size(x)));
n = length(x);
c = n/2+1;
y(2:(c-1)) = sqrt(2) .* real(y(2:(c-1)));
y((c+1):n) = -sqrt(2) .* imag(y((c+1):n));%

%
% Copyright (c) 2006. Iddo Drori
%

% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
