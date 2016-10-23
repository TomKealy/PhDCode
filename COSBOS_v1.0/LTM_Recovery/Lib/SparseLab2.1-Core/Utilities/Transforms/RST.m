function S = RST2(X)
% RST2: Real Sinusoid Transform of an n-vector
%  Usage:
%    S = RST(X);
%  Inputs:
%    X      input n-vector 
%  Outputs:
%    S      output vector which contains the transform coeffs. 
%
% Description
%   RST computes a 1-D real sinusoid transform of an n vector, n dyadic,
%   by taking the fft of X and using conjugate symmetry to eliminate
%   complex values. 
% See Also
%   Inv_RST, RST2

S = fft_mid0(X) ./ sqrt(prod(size(X)));

n = length(X);
n2 = n/2 + 1;    % Center point

S(2:(n2-1)) = sqrt(2) .* real(S(2:(n2-1)));
S((n2+1):n) = -sqrt(2) .* imag(S((n2+1):n));

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
