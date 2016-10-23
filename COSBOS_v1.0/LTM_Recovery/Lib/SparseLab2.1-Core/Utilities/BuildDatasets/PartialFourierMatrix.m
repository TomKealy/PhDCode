function A = PartialFourierMatrix(n,p)
% PartialFourierMatrix: creates random fourier matrix
% Usage
%	A = PartialFourierMatrix(n,p)
% Input
%	n           number of columns
%	p           number of rows
%                       
% Outputs
%	 A          random fourier matrix
%
F = fft(eye(p));
r = randperm(p);
A = F(r(1:n),:);
for j = 1:p
    A(:,j) = A(:,j) ./ norm(A(:,j));
end
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
