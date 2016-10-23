function X = Inv_RST2(S)
% Inv_RST2: Inverse Real Sinusoid Transform of an nxn array
%  Usage:
%    X = Inv_RST2(S);
%  Inputs:
%    S      n by n image 
%  Outputs:
%    X      reconstructed array. 
%
% Description
%   Inv_RST2 computes the inverse of RST2 by rearranging frequencies and 
%   taking the inverse 2-D fft.
% See Also
%   RST2

n = size(S,1);
n2 = n/2 + 1;    % Center point
isqrt2 = 1./sqrt(2);

X = S;

% Top row
X(1,2:(n2-1)) = isqrt2 .* (S(1,2:(n2-1)) + i .* S(1,n:-1:(n2+1)));
X(1,n:-1:(n2+1)) = isqrt2 .* (S(1,2:(n2-1)) - i .* S(1,n:-1:(n2+1)));

% Left Column
X(2:(n2-1),1) = isqrt2 .* (S(2:(n2-1),1) + i .* S(n:-1:(n2+1),1));
X(n:-1:(n2+1),1) = isqrt2 .* (S(2:(n2-1),1) - i .* S(n:-1:(n2+1),1));

% Middle column
X(2:(n2-1),n2) = isqrt2 .* (S(2:(n2-1),n2) + i .* S(n:-1:(n2+1),n2));
X(n:-1:(n2+1),n2) = isqrt2 .* (S(2:(n2-1),n2) - i .* S(n:-1:(n2+1),n2));

% Rest
X(2:n,2:(n2-1)) = isqrt2 .* (S(2:n,2:(n2-1)) + i .* S(n:-1:2,n:-1:(n2+1)));
X(n:-1:2,n:-1:(n2+1)) = isqrt2 .* (S(2:n,2:(n2-1)) - i .* S(n:-1:2,n:-1:(n2+1)));

X = real(ifft2_mid0(X)) .* sqrt(prod(size(X)));
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
