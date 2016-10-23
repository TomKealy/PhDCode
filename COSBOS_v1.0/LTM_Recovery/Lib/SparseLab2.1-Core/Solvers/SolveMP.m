function [sols, iters, activationHist] = SolveMP(A, y, N, maxIters, lambdaStop, solFreq, verbose, OptTol)
% SolveMP: Matching Pursuit
% Usage
%	[sols, iters, activationHist] = SolveMP(A, y, N, maxIters, lambdaStop, solFreq, verbose, OptTol)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%	maxIters    maximum number of iterations to perform. If not
%               specified, runs to stopping condition (default)
%   lambdaStop  If specified, the algorithm stops when the last coefficient 
%               entered has residual correlation <= lambdaStop. 
%   solFreq     if =0 returns only the final solution, if >0, returns an 
%               array of solutions, one every solFreq iterations (default 0). 
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default)
%	OptTol      Error tolerance, default 1e-5
% Outputs
%	 sols            solution(s) of MP
%    iters           number of iterations performed
%    activationHist  Array of indices showing elements entering  
%                    the solution set
% Description
%   SolveMP is an implementation of the Matching Pursuit scheme to 
%   estimate the solution of the sparse approximation problem
%      min ||x||_0 s.t. A*x = b
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as an m-file. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
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
    solFreq = 0;
end
if nargin < 5,
    lambdaStop = 0;
end
if nargin < 4,
    maxIters = 100*length(y);
end

explicitA = ~(ischar(A) || isa(A, 'function_handle'));
n = length(y);

% Initialize
x = zeros(N,1);
k = 1;
activeSet = [];
sols = [];
res = y;
normy = norm(y);
resnorm = normy;             
done = 0;

while ~done
    if (explicitA)
        corr = A'*res;
    else
        corr = feval(A,2,n,N,res,1:N,N); % = A'*y
    end
    [maxcorr i] = max(abs(corr));
    newIndex = i(1);
    
    if (explicitA)
        newVec = A(:,newIndex);
    else
        e = zeros(N,1);
        e(newIndex) = 1;
        newVec = feval(A,1,n,N,e,1:N,N);
    end
    
    % Project residual onto maximal vector
	newCoef = res' * newVec;
	res = res - newCoef .* newVec;
    resnorm = norm(res);
    
    % Update solution
    x(newIndex) = x(newIndex) + newCoef;
    activeSet = [activeSet newIndex];

    if ((resnorm <= OptTol*normy) | ((lambdaStop > 0) & (maxcorr <= lambdaStop)))
        done = 1;
    end

    if verbose
        fprintf('Iteration %d: Adding variable %d\n', k, newIndex);
    end

    k = k+1;
    if k >= maxIters
        done = 1;
    end

    if done | ((solFreq > 0) & (~mod(k,solFreq)))
        sols = [sols x];
    end
end

iters = k;
activationHist = activeSet;

%
% Copyright (c) 2006. Iddo Drori and Yaakov Tsaig
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
