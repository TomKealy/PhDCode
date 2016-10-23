function S = RST2(X)
% RST2: Real Sinusoid Transform of an nxn array
%  Usage:
%    S = RST2(X);
%  Inputs:
%    X      n by n image 
%  Outputs:
%    S      array which contains the transform coeffs. 
%
% Description
%   RST2 computes a 2-D real sinusoid transform of an nxn array, n dyadic,
%   by taking the 2-D fft of X and using conjugate symmetry to eliminate
%   complex values. 
% See Also
%   Inv_RST2

S = fft2_mid0(X) ./ sqrt(prod(size(X)));

n = size(X,1);
n2 = n/2 + 1;    % Center point
sqrt2 = sqrt(2);

% Top row
S(1,2:(n2-1)) = sqrt2 .* real(S(1,2:(n2-1)));
S(1,(n2+1):n) = -sqrt2 .* imag(S(1,(n2+1):n));

% Left Column
S(2:(n2-1),1) = sqrt2 .* real(S(2:(n2-1),1));
S((n2+1):n,1) = -sqrt2 .* imag(S((n2+1):n,1));

% Middle column
S(2:(n2-1),n2) = sqrt2 .* real(S(2:(n2-1),n2));
S((n2+1):n,n2) = -sqrt2 .* imag(S((n2+1):n,n2));

% Rest
S(2:n,2:(n2-1)) = sqrt2 .* real(S(2:n,2:(n2-1)));
S(2:n,(n2+1):n) = -sqrt2 .* imag(S(2:n,(n2+1):n));




%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
