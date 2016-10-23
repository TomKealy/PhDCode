function X = Inv_RST(S)
% Inv_RST: Inverse Real Sinusoid Transform of an n vector
%  Usage:
%    X = Inv_RST(S);
%  Inputs:
%    S      n vector 
%  Outputs:
%    X      reconstructed vector. 
%
% Description
%   Inv_RST computes the inverse of RST by rearranging frequencies and 
%   taking the inverse 1-D fft.
% See Also
%   RST

n = length(S);
n2 = n/2 + 1;    % Center point
isqrt2 = 1./sqrt(2);

X = S;

% Top row
X(2:(n2-1)) = isqrt2 .* (S(2:(n2-1)) + i .* S(n:-1:(n2+1)));
X(n:-1:(n2+1)) = isqrt2 .* (S(2:(n2-1)) - i .* S(n:-1:(n2+1)));

X = real(ifft_mid0(X)) .* sqrt(prod(size(X)));
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
