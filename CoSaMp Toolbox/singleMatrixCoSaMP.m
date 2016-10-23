function [alpha supp iters] = singleMatrixCoSaMP(y,Phi,k,opts)
%% Given measurements y = Phi*alpha, recover alpha via CoSaMP
% 
% Usage: [alpha supp iters] = singleMatrixCoSaMP(y,Phi,k,opts)
% Inputs: y - measurements of alpha
%         Phi - measurement matrix
%         k - number of nonzeros to estimate
%         opts - argument for setting various parameters
%           opts.maxiter - maximum number of iterations for CoSaMP
%           opts.tol - stopping criterion for CoSaMP 
%           opts.normBound - upper bound on norm of alpha
%                            (used for Tikhonov regularized least-squares)
%           opts.verbose - flag to show iteration progress
%
% Outputs: alpha - recovered signal
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
n = size(Phi,2);

if(isfield(opts,'maxiter'))
    maxiter = opts.maxiter;
else
    maxiter = 25;
end

if(isfield(opts,'tol'))
    tol = opts.tol;
else
    tol = 1e-6;
end

if(isfield(opts,'normBound'))
    normBound = opts.normBound;
else
    normBound = 2*norm(y);
end

if(isfield(opts,'verbose'))
    verbose = opts.verbose;
else
    verbose = 0;
end

%% Initialize CoSaMP variables
residual = y;
y_norm = norm(y);
residual_norm = 1;
iters = 0;
T = [];
supp = [];
TDiff = 1;
kMax = length(y);

%% Main loop of CoSaMP
while (iters <= maxiter) && (residual_norm>tol) && TDiff
    % Update iteration counter, support holder
    iters = iters + 1;
    oldT = T;     
    
    % Compute proxy = Phi^* residual
    proxy = ctranspose(Phi)*residual; 
    
    % Estimate support of proxy
    [junk,ind] = sort(abs(proxy),'descend');
    omega = ind(1:2*k);
    
    % Union support of proxy (trimmed if too large) with previously
    % estimated support
    [junk, ind] = setdiff(omega,supp); 
    omegaTrim = omega(ind); 
    T = union(omegaTrim(1:min(kMax-size(supp),end)),supp); 
    
    % Form estimate of alpha from y and submatrix corresponding to current
    % estimate of T
    PhiSub = Phi(:,T);
    alphaHat = lsOnSupport(y,PhiSub,normBound); 
    
    % Prune support estimate
    [junk,ind] = sort(abs(alphaHat),'descend');
    supp = T(ind(1:k));
    alphaHatTrim = alphaHat(ind(1:k));
        
    % Update residual and stopping criteria
    residual = y - Phi(:,supp)*alphaHatTrim; 
    residual_norm = norm(residual)/y_norm; 
    TDiff = length(union(T,oldT)) - length(intersect(T,oldT)); 
    
    % Display progress
    if(verbose)
        disp(['Iteration: ' num2str(iters) ' Norm(alphaHat): ' num2str(norm(alphaHat)) '   Residual: ' num2str(residual_norm) '  supp: ' mat2str(sort(supp))]);
    end
end

% Compute n-dimensional vector alpha
alpha = zeros(n,1);
alpha(supp) = alphaHatTrim;