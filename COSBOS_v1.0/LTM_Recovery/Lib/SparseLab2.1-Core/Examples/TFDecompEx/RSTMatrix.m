function Y = RSTMatrix(n)

I = eye(n);
for j = 1:n
    Y(:,j) = RST(I(:,j));
end

%
% Copyright (c) 2006. Iddo Drori
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
