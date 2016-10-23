function [x supp iters] = dpssCosamp(y,senseOpts,K,B,opts)
%% Given measurements y of a multiband signal x, recover x 
% via a modified version of CoSaMP
% 
% Usage: [x alpha supp iters] = dpssCosamp(y,senseOpts,K,B,opts)
% Inputs: y - measurements of x
%         senseOpts - options specifying sensing architecture
%         K - number of active bands
%         B - bandwidth of each band 
%             (B must be <= .5 and 1/B must be an integer)
%         opts - argument for setting various parameters
%           opts.N - dimension of x
%           opts.maxiter - maximum number of iterations for CoSaMP
%           opts.tol - stopping criterion for CoSaMP 
%           opts.dpssWidth - width of DPSS basis vectors to use in recovery
%           opts.dpssSize - number of DPSS basis vectors to use in recovery
%           opts.PsiRecover - precomputed DPSS basis for faster running
%           opts.PsiDetect - precomputed DPSS basis for faster running
%           opts.modulator - precomputed modulator matrix for faster running
%           opts.normBound - upper bound on norm of x (used for Tikhonov
%                            regularize least-squares)
%           opts.detectAlg - algorithm to use in computing projection onto
%                            the set of block-sparse signals
%           opts.verbose - flag to show iteration progress
%
% Outputs: x - recovered signal
%          supp - support of recovered signal
%          iters - number of iterations required to converge
%
% Most recent change - 9/8/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.


%% Set optional parameters
if(isfield(opts,'N'))
    N = opts.N;
else
    N = length(sense(y,1,senseOpts));
end

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

if(isfield(opts,'dpssWidth'))
    dpssWidth = opts.dpssWidth;
else
    dpssWidth = B*N/2;
end

if(isfield(opts,'dpssSize'))
    dpssSize = opts.dpssSize;
else
    dpssSize = 2*dpssWidth;
end

if(isfield(opts,'PsiRecover'))
    PsiRecover = opts.PsiRecover;
else
    PsiRecover = dpss(N,dpssWidth,dpssSize); 
end

if(isfield(opts,'PsiDetect'))
    PsiDetect = opts.PsiDetect;
else
    PsiDetect = PsiRecover; 
end

if(isfield(opts,'modulator'))
    modulator = opts.modulator;
else
    modulator = generateModulator(N,B); 
end

if(isfield(opts,'normBound'))
    normBound = opts.normBound;
else
    normBound = 2*norm(y);
end

if(isfield(opts,'detectAlg'))
    detectOpts.alg = opts.detectAlg;
else
    detectOpts.alg = 'blockOMP'; 
end

if(isfield(opts,'debias'))
    debias = 1;
else
    debias = 0;
end

if(isfield(opts,'verbose'))
    verbose = opts.verbose;
else
    verbose = 0;
end


%% Initialize CoSaMP variables
x = zeros(N,1);
residual = y;
y_norm = norm(y);
residual_norm = 1;
iters = 0;
T = [];
supp = [];
TDiff = 1;
KMax = floor(length(y)/size(PsiRecover,2));

%% Main loop of CoSaMP
while (iters <= maxiter) && (residual_norm>tol) && TDiff
    % Update iteration counter, support holder
    iters = iters + 1;
    oldT = T;     
    
    % Compute proxy = Phi^T residual
    proxy = sense(residual,1,senseOpts); 
    
    % Estimate support of proxy
    omega = bandDetect(proxy,2*K,PsiDetect,PsiRecover,modulator,detectOpts); 
    
    % Union support of proxy (trimmed if too large) with previously
    % estimated support
    [junk, ind] = setdiff(omega,supp); 
    omegaTrim = omega(ind); 
    T = union(omegaTrim(1:min(KMax-size(supp),end)),supp); 
    
    % Form signal estimate from y and currently estimated support
    xHat = bandRecover(y,senseOpts,T,PsiRecover,modulator,normBound); 
    
    % Estimate support of xHat and compute pruned x
    [supp, x] = bandDetect(xHat,K,PsiDetect,PsiRecover,modulator,detectOpts); 
    
    if(debias)
        x = bandRecover(y,senseOpts,supp,PsiRecover,modulator,normBound);
    end
           
    % Update residual and stopping criteria
    residual = y - sense(x,0,senseOpts); 
    residual_norm = norm(residual)/y_norm; 
    TDiff = length(union(T,oldT)) - length(intersect(T,oldT)); 
    
    % Display progress
    if(verbose)
        disp(['Iteration: ' num2str(iters) '   Residual: ' num2str(residual_norm) '  T: ' mat2str(sort(T))]);
    end
end