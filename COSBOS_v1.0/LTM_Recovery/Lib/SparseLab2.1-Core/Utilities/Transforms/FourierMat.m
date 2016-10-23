function F = FourierMat(N);
% This function generates an NxN Fourier Transform matrix F, s.t. writing
% x_hat = F*x is equivalent to x_hat = fft(x) ./ sqrt(N)

[n, k] = meshgrid(0:(N-1));
F = exp(-sqrt(-1) .* 2.*pi .* k .* n ./ N);
F = F ./ sqrt(N); 
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
