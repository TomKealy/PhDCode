function I = GenRandomRects(n,k);
% This function generates an nxn image of k random rectangles

I = zeros(n);

for j = 1:k
    x1 = floor(rand * n);
    x2 = floor(x1 + rand * (n - x1));
    y1 = floor(rand * n);
    y2 = floor(y1 + rand * (n - y1));
    I(y1:y2,x1:x2) = I(y1:y2,x1:x2) + rand * ones(y2-y1+1, x2-x1+1);
end

    
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
