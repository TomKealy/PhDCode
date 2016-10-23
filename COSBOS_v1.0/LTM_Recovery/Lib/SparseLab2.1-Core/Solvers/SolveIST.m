function [sols, numIters, activeSet] = SolveIST(A, y, p, mu, gamma, maxIters, fullPath, verbose, OptTol)
% SolveIST: Iterative Soft Thresholding by threshold path following
% Usage
%	[sols, numIters, activationHist] = SolveIST(A, y, p, mu, gamma, maxIters, fullPath, verbose, OptTol)
% Input
%	A           nxp matrix
%	y           n-vector
%   p           size(A,2)
%   gamma       step size
%   mu          threshold update, default 1-gamma
%	maxIters    
%   fullPath    1 returns entire solution path, 0 final
%   verbose     
%	OptTol      error tolerance, default 1e-5
% Outputs
%	 sols       solution of BPDN problem
%	 numIters  
% Description
%   SolveIST executes the Iterative Soft Thresholding scheme by least squares update,
%   approximates LARS.
%
fastop = (ischar(A) || isa(A, 'function_handle'));

n = length(y);

if nargin < 9,
    OptTol = 1e-5;
end
if nargin < 8,
    verbose = 1;
end
if nargin < 7,
    fullPath = 0;
end
if nargin < 6,
    maxIters = 50;
end
if nargin < 5,
    gamma = 1/2;
end
if nargin < 4,
    mu = 1-gamma;
end
if nargin < 3,
    p = size(A,2);
end

x = zeros(p,1);
activeSet = [];
iter = 1;
sols = [];
epsilon = 1e-6;

res = y;
if fastop
    corr = feval(A,2,p,res);
else
    corr = A'*res;
end

normy = norm(y);
lambda = max(abs(corr));

% lsqrms params
damp = 0;
atol = 1e-6;
btol = 1e-6;
conlim = 1e+10;
itnlim = 10;
show = 0;

while  (iter <= maxIters) & (norm(res) > OptTol*normy)

    activeSet = find(abs(corr) > lambda - epsilon)';

    if (verbose)
        fprintf('Iteration %d: |I| = %d, lambda = %g, ||r|| = %g\n', iter, length(activeSet), lambda, norm(res));
    end

    dx = zeros(p,1);
    if fastop
        [dx(activeSet), istop, itn] = lsqrms(n, length(activeSet), 'FastOpLS', activeSet, p, res, damp, atol, btol, conlim, itnlim, show);
    else
        AI = A(:,activeSet);
        dx(activeSet) = AI \ res;
    end

    nres = res;
    gamma = 1;
    xt = x;
    while norm(nres) >= norm(res) & gamma > 0.1
        x = xt + gamma*dx;
        nres = y - feval(A,1,p,x);
        gamma = gamma-0.1;
    end
        
    if fastop
        res = y - feval(A,1,p,x);
        corr = feval(A,2,p,res);
    else
        res = y - A*x;
        corr = A'*res;
    end

    lambda = mu*lambda;

    if fullPath
        sols = [sols x];
    end

    iter = iter+1;
end

if ~fullPath
    sols = x;
end

numIters = iter;

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
