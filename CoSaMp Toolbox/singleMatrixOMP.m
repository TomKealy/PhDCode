function [alpha supp] = singleMatrixOMP(y,Phi,k,opts)
%% Given measurements y = Phi*alpha, recover alpha via OMP
% 
% Usage: [alpha supp] = singleMatrixOMP(y,Phi,k,opts)
% Inputs: y - measurements of alpha
%         Phi - measurement matrix  
%         k - number of nonzeros to estimate
%         opts - argument for setting various parameters
%           opts.normBound - upper bound on norm of alpha 
%           (used for Tikhonov regularized least-squares)
%           opts.verbose - flag to show iteration progress
%
% Outputs: alpha - recovered signal
%          supp - support of recovered signal
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

%% Initialize OMP variables
residual = y;
supp = [];

%% Main loop of OMP
for jj=1:k,
    
    % Compute proxy = Phi^* residual
    proxy = ctranspose(Phi)*residual; 
    
    % Find the elemnt of proxy with largest magnitude
    energy = abs(proxy);
    energy(supp) = 0;
    [junk,ind] = sort(energy,'descend');
    supp = union(supp,ind(1));
    
    % Form estimate of alpha from y and submatrix corresponding to current
    % estimate of supp
    PhiSub = Phi(:,supp);
    alphaHat = lsOnSupport(y,PhiSub,normBound);     
    
    % Update residual and stopping criteria
    residual = y - PhiSub*alphaHat; 
    
    % Display progress
    if(verbose)
        disp(['Iteration: ' num2str(jj) ' Norm(alphaHat): ' num2str(norm(alphaHat)) '  supp: ' mat2str(sort(supp))]);
    end
end

% Compute n-dimensional vector alpha
alpha = zeros(n,1);
alpha(supp) = alphaHat;