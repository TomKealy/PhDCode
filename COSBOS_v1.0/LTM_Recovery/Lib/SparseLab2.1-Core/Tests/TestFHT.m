% TestFHT

% FHT.m is located in Utilities/Transforms/ 

n = 256;
A = eye(n);
B = zeros(n);
for ii = 1:n
    B(:,ii) = FHT(A(:,ii));
end

twonorm(B - hadamard(n) ./ sqrt(n))

x = rand(n,1);
y1 = FHT(x);
y2 = (hadamard(n) ./ sqrt(n)) * x;
twonorm(y1-y2)

x = rand(n,1);
y = rand(n,1);
Ax = FHT(x);
ATy = FHT(y);
twonorm(y'*Ax - x'*ATy)

x = rand(n,1);
y = FHT(x);
x2 = FHT(y);
twonorm(x-x2)

%
% Copyright (c) 2006. Yaakov Tsaig
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
