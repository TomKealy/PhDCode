function A = PartialHadamardMatrix(n,p)
% PartialHadamardMatrix: creates random Hadamard matrix
% Usage
%	A = PartialHadamardMatrix(n,p)
% Input
%	n           number of columns
%	p           number of rows
%                       
% Outputs
%	 A          random Hadamard matrix
%
H = hadamard(p);
r = randperm(p);
A = H(r(1:n),:);
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
