function [sol, numIters] = SolveStOMP(A, y, N, thresh, param, maxIters, verbose, OptTol)
% SolveStOMP: Implementation of Iterative Threshold-Selective Projection
% algorithm
% Usage
%	[sol, numIters] = SolveStOMP(A, y, N, thresh, param, maxIters,
%	verbose, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%   thresh      Thresholding strategy: FDR or FAR. default is FDR.
%   param       Sensitivity parameter for threshold selection.
%	maxIters    maximum number of StOMP iterations to perform, default 10. 
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default).
%	OptTol      Error tolerance, default 1e-5.
% Outputs
%	 sol        Solution of StOMP.
%	 numIters   Total number of steps taken.
% Description
%   SolveStOMP implements the Stagewise Ortogonal Matching Pursuit, 
%   as described in the paper
%   D.L. Donoho et al., "Sparse Solution of Underdetermined Linear 
%   Equations by Stagewise Ortogonal Matching Pursuit". 
%   It essentially computes an approximate solution to the problem 
%     min ||x||_1  s.t.  Ax = y
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as an m-file. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
%   This function is used as input to Michael Saunders' least-squares
%   solver, LSQR, available at:
%     http://www.stanford.edu/group/SOL/software/lsqr.html
% See Also
%   SolveLasso, SolveBP, SolveOMP
%

if nargin < 8,
    OptTol = 1e-5;
end
if nargin < 7,
    verbose = 0;
end
if nargin < 6,
    maxIters = 10;
end
if nargin < 4,
    thresh = 'FDR';
    param = 0.5;
end

% Initialize threshold parameters
switch lower(thresh)
    case 'fdr'
        q = param;
    case 'far'
        alpha_0 = param;
end

explicitA = ~(ischar(A) || isa(A, 'function_handle'));
n = length(y);

% Set LSQR parameters
damp   = 0;
atol   = 1.0e-6;
btol   = 1.0e-6;
conlim = 1.0e+10;
itnlim = n;
show   = 0;

% Initialize
iter = 1;
res = y;  
normy = norm(y);
activeSet = [];
x_I = [];
Ifull = 1:N;
done = 0;

while ~done
    % Compute residual correlations
    if (explicitA)
        corr = sqrt(n) .* A'*res ./ norm(res);
    else
        corr = feval(A,2,n,N,res,Ifull,N);
        corr = sqrt(n) .* corr ./ norm(res);
    end
    
    switch lower(thresh)
        case 'fdr'
            thr = fdrthresh(corr, q);
        case 'far'
            thr = norminv(1 - alpha_0/2, 0, 1);
    end

    % Apply hard thresholding
    I = find(abs(corr) > thr);
    
    % If no new entries, we are done
    J = union(activeSet, I);
    if (length(J) == length(activeSet)) 
        done = 1;
    else
        activeSet = J;
        
        % Compute current estimate and residual
        if (explicitA)
            x_I = A(:, activeSet) \ y;
            res = y - A(:, activeSet)*x_I;
        else
            % Solve (A_I^T*A_I)x_I = y using lsqr
            p = length(activeSet);
            %[x_I, istop, itn] = lsqr_ms(n, p, A, activeSet, N, y, damp, atol, btol, conlim, itnlim, show);
            [x_I, flag] = lsqr(@lsqrMat,y,OptTol,20);

            % Compute residual
            Ax_I = feval(A,1,n,p,x_I,activeSet,N);
            res = y - Ax_I;
        end

        if (verbose)
            fprintf('Iteration %d: |I| = %d, ||r||_2 = %g\n', iter, length(activeSet), norm(res));
        end
    end
    
    iter = iter+1;
    
    % Check stopping criteria
    if (iter > maxIters)
        done = 1;
    end
    if norm(res) <= OptTol*normy
        done = 1;
    end
end

sol = zeros(N,1);
sol(activeSet) = x_I;
numIters = iter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function y = lsqrMat(x,transpose)
        switch transpose
            case 'transp'
                y = feval(A,2,n,p,x,activeSet,N);
            case 'notransp'
                y = feval(A,1,n,p,x,activeSet,N);
        end
    end

end

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
