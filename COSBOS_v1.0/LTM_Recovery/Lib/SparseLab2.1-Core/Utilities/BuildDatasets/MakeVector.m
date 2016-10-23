function x = MakeVector(k, p, type)
% MakeVector: creates various types of sparse test vectors
% Usage
%	x = MakeVector(k, p, type)
% Input
%	k           sparsity level
%	p           vector length
%	type        one of: uniform, gaussian, cauchy, sign
%                       
% Outputs
%	 x          1xp test vector
%
if nargin < 3
    type = 'uniform';
end

x = zeros(p,1);
if strcmp(type, 'uniform')
    x(1:k) = rand(k,1);
elseif strcmp(type, 'gaussian')
    x(1:k) = randn(k,1);
elseif strcmp(type, 'cauchy')
    x(1:k) = randn(k,1)./randn(k,1);
elseif strcmp(type, 'sign')
    x(1:k) = sign(randn(k,1));
end


%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
