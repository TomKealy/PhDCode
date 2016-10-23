function [x alpha supp iters] = doubleMatrixCoSaMP(y,A,D,k,opts)
%% Given measurements y of a signal x, recover x via Signal Space CoSaMP
%
% Usage: [x alpha supp iters] = doubleMatrixCoSaMP(y,A,D,k,opts)
% Inputs: y - measurements of x
%         A - sensing matrix A
%         D - sparsity dictionary D
%           note that A and D are required as separate inputs; hence the
%           filename "doubleMatrixCoSaMP"
%         k - number of nonzeros to estimate
%         opts - argument for setting various parameters
%           opts.maxiter - maximum number of iterations for CoSaMP
%           opts.tol - stopping criterion for CoSaMP
%           opts.normBound - upper bound on norm of alpha
%                            (used for Tikhonov regularized least-squares)
%           opts.verbose - flag to show iteration progress
%
% Outputs: x - recovered signal
%          alpha - recovered coefficient vector
%          supp - support of recovered signal
%          iters - number of iterations required to converge
%
% Most recent change - 7/31/2012
%
% Copyright 2012, Mark Davenport, Deanna Needell, Michael Wakin
%
% This file is part of Signal Space CoSaMP Toolbox version 1.0.
%
%    Signal Space CoSaMP Toolbox is free software: you can redistribute it 
%    and/or modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation, either version 3 of the 
%    License, or (at your option) any later version.
%
%    Signal Space CoSaMP Toolbox is distributed in the hope that it will be 
%    useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with Signal Space CoSaMP Toolbox.  If not, 
%    see <http://www.gnu.org/licenses/>.

%% Set parameters
[n d] = size(D);

if(isfield(opts,'maxiter'))
    maxiter = opts.maxiter;
else
    maxiter = 50;
end

if(isfield(opts,'tol'))
    tol = opts.tol;
else
    tol = 1e-6;
end

if(isfield(opts,'verbose'))
    verbose = opts.verbose;
else
    verbose = 0;
end

if(isfield(opts,'normBound'))
    normBound = opts.normBound;
else
    normBound = 1e10;
end
singleMatrixOpts.normBound = normBound;

%% Initialize CoSaMP variables
residual = y;
y_norm = norm(y);
residual_norm = 1;
iters = 0;
T = [];
supp = [];
TDiff = 1;
kMax = length(y);

if strcmp(opts.alg,'nnd')
    Dnorms = sqrt(diag(D'*D));
end

%% Main loop of CoSaMP
while (iters <= maxiter) && (residual_norm>tol) && TDiff
    % Update iteration counter, support holder
    iters = iters + 1;
    oldT = T;
    
    % Compute proxy = A^* residual
    proxy = ctranspose(A)*residual;
    
    % Estimate support of proxy
    if strcmp(opts.alg,'omp')
        [junk,omega] = singleMatrixOMP(proxy,D,2*k,singleMatrixOpts);
    elseif strcmp(opts.alg,'cosamp')
        singleMatrixOpts.maxiter = 2*maxiter;
        [junk,omega] = singleMatrixCoSaMP(proxy,D,2*k,singleMatrixOpts);
    elseif strcmp(opts.alg,'l1')
        spgOpts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
        alphaHat = spg_bp(D, proxy, spgOpts);         % Solves alphaHat = argmin ||alpha||_1 s.t. D*alpha = proxy
        [alphaSort,alphaInds] = sort(abs(alphaHat),'descend');
        omega = alphaInds(1:2*k);
    elseif strcmp(opts.alg,'nnd')
        % Use this only if D is a non-normalized orthogonal dictionary
        alphaProxy = ctranspose(D)*proxy./Dnorms;
        [junk,ind] = sort(abs(alphaProxy),'descend');
        omega = ind(1:2*k);
    end
    
    % Union support of proxy (trimmed if too large) with previously
    % estimated support
    [junk, ind] = setdiff(omega,supp);
    omegaTrim = omega(ind);
    T = union(omegaTrim(1:min(kMax-size(supp),end)),supp);
    
    % Form signal estimate from y and currently estimated support
    basis = D(:,T);
    projBasis = A*basis;
    alphaHat = lsOnSupport(y,projBasis,normBound);
    xHat = basis*alphaHat;
    
    % Estimate support of xHat and compute pruned x
    if strcmp(opts.alg,'omp')
        [alpha,supp] = singleMatrixOMP(xHat,D,k,singleMatrixOpts);
    elseif strcmp(opts.alg,'cosamp')
        singleMatrixOpts.maxiter = maxiter;
        [alpha,supp] = singleMatrixCoSaMP(xHat,D,k,singleMatrixOpts);
    elseif strcmp(opts.alg,'l1')
        spgOpts = spgSetParms('verbosity',0);         % Turn off the SPGL1 log output
        alphaHat = spg_bp(D, xHat, spgOpts);          % Solves alphaHat = argmin ||alpha||_1 s.t. D*alpha = xHat
        [alphaSort,alphaInds] = sort(abs(alphaHat),'descend');
        supp = alphaInds(1:k);        
        alpha = zeros(d,1);
        alpha(supp) = lsOnSupport(xHat,D(:,supp),normBound); % Debiasing step after L1
    elseif strcmp(opts.alg,'nnd')
        % Use this only if D is a non-normalized orthogonal dictionary
        alphaProxy = ctranspose(D)*xHat./Dnorms;
        [junk,ind] = sort(abs(alphaProxy),'descend');
        supp = ind(1:k);
        alpha = zeros(d,1);
        alpha(supp) = lsOnSupport(xHat,D(:,supp),normBound);
    end
    x = D*alpha;
    
    % Update residual and stopping criteria
    residual = y - A*x;
    residual_norm = norm(residual)/y_norm;
    TDiff = length(union(T,oldT)) - length(intersect(T,oldT));
    
    % Display progress
    if(verbose)
        disp(['Iteration: ' num2str(iters) ' Norm(x): ' num2str(norm(x)) ' Norm(alpha): ' num2str(norm(alpha)) ' Residual: ' num2str(residual_norm)]);
        disp(['  Omega: ' mat2str(sort(omega))]);
        disp(['  T: ' mat2str(sort(T))]);
        disp(['  supp: ' mat2str(sort(supp))]);
    end
end
