function y = Inv_FHT2(x)
% Inv_FHT2: Computes the inverse 2-D Hadamard transform of a signal
%  Usage:
%    y = Inv_FHT2(x);
%  Inputs:
%    x      Hadamard transform coefficients
%  Outputs:
%    y      image
%
% Description
%   Inv_FHT 2applies the forward 2-D Hadamard transform, since it is a self-adjoint operator.
% See Also
%   FHT2

y = FHT2(x);
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
