function [x, cost] = AMP(A, y, tau, xstar, MAX_ITER)

t_start = tic;

[m, n] = size(A);
k = 1/m;

x = A'*y;
zt = y;

cost = zeros(1, MAX_ITER);

for i=1:MAX_ITER
    gammat = tau*k*norm(zt)^2;
    x = shrinkage(x + A'*zt, gammat);
    zt = y - A*x + k*nnz(x)*zt;
    cost(i) = norm(x - xstar)/norm(xstar);
end

toc(t_start);

end

function z = shrinkage(x, kappa)
    z = max( 0, x - kappa ) - max( 0, -x - kappa );
end