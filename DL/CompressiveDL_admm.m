function [xhat, D] = CompressiveDL_admm(B, A, D0, Y, sigma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

max_iter = 100;

%init D to be random

[m, n] = size(A);
D = D0;
p = norm(D0, 'fro');
lambda = sigma*sqrt(2*log(n))/100;

max_eig = max(abs(eig((A*D)'*(A*D))));

rho = nthroot(1/max_eig, 3);

[xhat] = lasso_dl(A*D, B, lambda, rho, 1.0);
i=0;
for i =1:100
    i
    [xhat] = lasso_dl(A*D, B, lambda, rho, 1.0);
    G = -A'*(B-A*D*xhat)*xhat';
    mu = norm(G, 'fro')/norm(A*G*xhat)
    D = D + mu*G;
    D = sqrt(p)*D/norm(D, 'fro');
end

end