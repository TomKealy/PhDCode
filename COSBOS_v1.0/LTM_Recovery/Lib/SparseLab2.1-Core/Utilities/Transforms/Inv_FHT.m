function y = Inv_FHT(x)
% Inv_FHT: Computes the inverse Hadamard transform of a signal
%  Usage:
%    y = Inv_FHT(x);
%  Inputs:
%    x      Hadamard transform coefficients
%  Outputs:
%    y      signal
%
% Description
%   Inv_FHT applies the forward Hadamard transform, since it is a self-adjoint operator.
% See Also
%   FHT

y = FHT(x);
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
