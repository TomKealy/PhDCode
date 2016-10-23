% Test the routine SolveStOMP

n = 400;
d = 200;

err = [];
for k = 5:5:100
    A = MatrixEnsemble(d,n);
    x = 100*[rand(k,1)-0.5*ones(k,1); zeros(n-k,1)];
    y = A*x;

    % Initialize threshold parameters
    delta = d/n;
    rho = k/d;
    S = 10;
    alpha_0 = delta*(1-rho)/S;
    q = min((d-k)/k,0.5);

    [sol, iters] = SolveStOMP(A, y, n, 'FAR', alpha_0, S, 1);
    err = [err norm(sol - x)];
end

plot(5:5:100,err);
    
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
