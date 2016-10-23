function [x supp] = fftCosamp(y,senseOpts,K,opts)
%% Given measurements y of a signal x, recover x 
% via CoSaMP and the DFT basis
% 
% Usage: [x supp] = fftCosamp(y,pnseq,K,opts)
% Inputs: y - measurements of x
%         senseOpts - options specifying sensing architecture
%         K - sparsity
%         opts - argument for setting various parameters
%           opts.maxiter - maximum number of iterations for CoSaMP
%           opts.tol - stopping criterion for CoSaMP 
%
% Outputs: x - recovered signal
%          supp - support of recovered signal
%
% Most recent change - 9/10/2011
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

%% Set signal parameters
N = opts.N;
M = length(y);

%% Set optional parameters
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

%% Initialize CoSaMP variables
fft_basis = sqrt(N)*ifft(eye(N));  % Generate explicit DFT matrix
x = zeros(N,1);
residual = y;
y_norm = norm(y);
residual_norm = 1;
iters = 0;
supp = [];

%% Main loop of CoSaMP
while (iters <= maxiter) && (residual_norm>tol)
    % Update iteration counter 
    iters = iters + 1; 
    
    % Compute proxy = Phi^T residual
    proxy = sense(residual,1,senseOpts); 
    
    % Calculate omega
    PROXY = 1/sqrt(N)*fft(proxy);  % Compute fft
    [junk,index] = sort(abs(PROXY),'descend');  % Sort fft coefficients
    omega = index(1:2*K);  % Find 2K largest elements
    
    % Union support of proxy with previously estimated support
    [junk, ind] = setdiff(omega,supp); 
    omegaTrim = omega(ind); 
    T = union(omegaTrim(1:min(M-size(supp),end)),supp);
    
    % Form signal estimate from y and currently estimated support
    projBasis = zeros(M,length(T));
    for jj=1:length(T),
       projBasis(:,jj) = sense(fft_basis(:,T(jj)),0,senseOpts); 
    end
    alpha = zeros(N,1);
    alpha(T) = projBasis\y;
    
    % Calculate pruned support
    [junk,index] = sort(abs(alpha),'descend');  % Sort fft coefficients
    supp = index(1:K);  % Find K largest elements
    
    % Compute pruned x
    alpha_prune = zeros(N,1);
    alpha_prune(supp) = alpha(supp);
    x = sqrt(N)*ifft(alpha_prune);
    
    % Update residual and stopping criterion
    residual = y - sense(x,0,senseOpts);     
    residual_norm = norm(residual)/y_norm; 
end