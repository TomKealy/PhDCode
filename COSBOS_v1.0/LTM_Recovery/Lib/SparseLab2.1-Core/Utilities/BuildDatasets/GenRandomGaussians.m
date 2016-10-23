function I = GenRandomGaussians(n,k,state);
% This function generates an nxn image of k random Gaussians, using state
% as the random seed.

rand('state',state);

I = zeros(n);
gSize = 30;
meanAmp = 10;

for j = 1:k
    x1 = floor(rand * (n-gSize))+1;
    y1 = floor(rand * (n-gSize))+1;
    sigma = gSize * rand / 3;
    I(x1:(x1+gSize-1),y1:(y1+gSize-1)) = I(x1:(x1+gSize-1),y1:(y1+gSize-1)) + ...
        meanAmp*exprnd(1).*fspecial('gaussian',gSize,sigma);
end

    
%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
