function [sols, resids, numIters] = SolveIRWLS(A, y, p, maxIters, fullPath, verbose, OptTol)
% SolveIRWLS: Implementation of Iterative Reweighted Least Squares
% Usage
%	[sols, resids, numIters] = SolveIRWLS(A, y, maxIters, fullPath, verbose, OptTol)
% Input
%	A           nxp matrix
%	y           nx1 observation vector
%	maxIters 
%   fullPath    1 to return the entire solution path, 0 to return solution
%   verbose     1 to print out, 0 otherwise
%	OptTol      Error tolerance, default 1e-5
% Outputs
%	 sols       solution path
%    resids     residuals
%	 numIters   Total number of steps taken
% Description
%   SolveIRWLS executes the Iterative Reweighted Least Squares scheme
%   and solves min ||y-Ax||_1
% See Also
%   SolveIST, pdco.m, pdcoSet.m
%

n = length(y);

if nargin < 4
    maxIters = 50; % one iteration is LS solution
end

if nargin < 5
    fullPath = 0;
end

if nargin < 6
    verbose = 1;
end

if nargin < 7
    OptTol = 1e-8;
end

fastop = (ischar(A) || isa(A, 'function_handle'));

% lsqrms params
damp = 0;
atol = 1e-6;
btol = 1e-6;
conlim = 1e+10;
itnlim = 10;
show = 0;

iter = 1;
x = zeros(p,1);

if (fastop)
    [x, istop, itn] = lsqrms(n, p, 'FastOpLS', [1:p], p, y, damp, atol, btol, conlim, itnlim, show);
    res = feval(A,1,p,x);
else
    x = A \ y;
    res = y - A*x;
end

epsk = 1;
global W;
W = diag(1./abs(res + epsk));
sols = [];
resids = [];

while (iter <= maxIters) & (epsk >= OptTol)

    if (fastop)
        [x, istop, itn] = lsqrms(n, p, 'FastOpLSW', [1:p], p, W*y, damp, atol, btol, conlim, itnlim, show);
    else
        x = W*A \ W*y;
    end
    
    if (verbose)
        fprintf('Iteration %d: epsilon = %g\n', iter, epsk);
    end
    
    if (fastop)
        res = y - feval(A, 1, p, x);
    else
        res = y - A*x;
    end
    
    epsk = 1./sqrt(iter).*mean(abs(res));
    W = diag(1./abs(res + epsk));
    
    if fullPath
        sols = [sols x];
        resids = [resids res];
    end

    iter = iter+1;
end

if ~fullPath
    sols = x;
    resids = res;
end

numIters = iter-1;
    
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
