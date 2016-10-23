function [xhat, D, A] = SparseDL_admm(Y, B, A, beta, sigma, T, D0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

max_iter = 10;

%init D to be random

[m, n] = size(A);
[k, l] = size(Y);
D = D0; %randn(n, n);
lambda = sigma*sqrt(2*log(n))/100;
xhat = randn(n,l);
max_eig = max(abs(eig((A*D)'*(A*D))));

rho = nthroot(1/max_eig, 3);

for ii =1:max_iter %norm(Y-D*xhat, 'fro')>1 || ii > max_iter
    ii
    %norm(Y-D*xhat, 'fro')
    [xhat] = lasso_dl(A*D, B, lambda, rho, 1.0);
    A = (B*(D*xhat)')*inv((D*xhat)*(D*xhat)' + 2*beta*speye(n));
    D = (Y*xhat')*inv(xhat*xhat'+ 2*beta*speye(n));
end

end

