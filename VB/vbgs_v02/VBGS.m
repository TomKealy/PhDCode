function [w, varargout] = VBGS(y, A, grouping, options)
% Variational Bayesian Group-Sparse Estimation
% Author (a.k.a. person to blame): S. Derin Babacan
% Created: December 5, 2011
% Last updated: June 15,2012 (added dimensionality reduction)
%
% Usage:
% [w, varargout] = VBGS(y, A, grouping, options)
% Solves for w in
% y = Aw + n where w is group sparse and n is 
% dense Gaussian noise
%
% If you're using this code, please acknowledge 
%    S. D. Babacan, S. Nakajima, M. N. Do, 
%    "Bayesian Group-Sparse Modeling and Variational Inference," 
%    submitted to IEEE Transactions on Signal Processing, 2012.
% 
% -----------------------------------------------------------------------
%                                INPUTS
%   y             : Mx1 observation vector
%   A             : MxN dictionary (or measurement) matrix
%   grouping      : cell array of length G showing the group structure
%                   each element in grouping (grouping{i}) contains an array of 
%                   indices that belong to group i
% 
%               -------------------------------------------
%   options       : All options are *optional*. The default values are set
%                   automatically.
%                    
%   dist          : signal distribution. 
%                   'McKay', 'MLaplace', 'MStudent', 'Jeffreys' (default)
%
%   pt_estimate   : the type of point estimates 
%                   'Mode' or 'Mean' (default)
%
%   cov_met       : Covariance estimation
%                   1: Full covariance with O(N^3) estimation
%                   2: Full covariance with O(M^3) estimation (default)
%                   3: Approximated covariance with O(N) estimation
%                      (suggested in large scale problems)                    
%                   4: MAP: covariance = 0
% 
%   verbose       : output the progress? (0/1) default: 0. 
%   init          : Initialization method. 
%                   'rand': initialize w by a random vector
%                   'ml'  : initialize w by the maximum likelihood solution
%                           A\y  (default)
%   k_a
%   theta_a       : hyperparameter values for a. default: 1e-6
%
%   k_b
%   theta_b       : hyperparameter values for b. default: 1e-6
%
%   a_beta0
%   theta_beta0   : hyperparameter values for beta. default: 0
%   
%   beta_init     : initial value for noise precision beta
%                   either a double value or 'auto' (default)
%   
%   conv_thr      : convergence threshold (default: 1e-7)
%   MAXITER       : max number of iterations. default: 100
%   UPDATE_BETA   : update noise variance? (0/1). default: 1
%   w_true        : Original w vector for simulations. 
%                   If supplied and verbose=1,
%                   the error is reported at each iteration
% -----------------------------------------------------------------------
%                                OUTPUTS
% w               : estimated group sparse vector
% varargout{1}    : number of iterations
% varargout{2}    : estimated inverse variances z
% varargout{3}    : estiamted noise precision beta
% varargout{4}    : estimated covariance matrix Sigma_w
% varargout{5}    : hyperparameter a
% varargout{6}    : hyperparameter b
% -----------------------------------------------------------------------
%
% This code is not optimized for speed, but it should run fairly fast for
% small to medium scale problems. 
% Important: in this code, z denotes the INVERSE variances, while z are the
% variances in the paper.
%
% Copyright (c): S. Derin Babacan, 2011
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.
%


%% Check input variables, 
if ~exist('options','var')
    options = populate_vars([]);
else
    options = populate_vars(options);
end

G = length(grouping);

verbose                = options.verbose;
k_a                    = options.k_a;
k_b                    = options.k_b;
theta_a                = options.theta_a;
theta_b                = options.theta_b;
lambda                 = options.lambda;
MAXITER                = options.MAXITER;
UPDATE_BETA            = options.UPDATE_BETA;
UPDATE_BETA_init       = options.UPDATE_BETA_init;
k_beta0                = options.k_beta0;
theta_beta0            = options.theta_beta0;
ESTIMATE_HYPERPARAMS   = options.ESTIMATE_HYPERPARAMS;
dist                   = options.dist;
pt_estimate            = options.pt_estimate;
cov_met                = options.cov_met;
conv_thr               = options.conv_thr;
a0                     = options.a0;
b0                     = options.b0;

a                      = ones(G,1)*a0;
b                      = ones(G,1)*b0;
DIMRED                 = options.DIMRED;


[m n] = size(y);
[M N] = size(A);
if G == N,
    INDIVIDUAL_GROUPS = 1;
else
    INDIVIDUAL_GROUPS = 0;
end

Neff = N;

used_indices = [1:N];

if DIMRED,
    if isfield(options, 'DIMRED_THR')
        DIMRED_THR = options.DIMRED_THR;
    else
        DIMRED_THR = 1e4;
    end
end

%% ------------------------------------------------------------------------------
% Dimensionality reduction currently does NOT work with overlapping groups
% Consistency check
Nc = 0;
for i=1:G,
    Nc = Nc + length(grouping{i});
end
if Nc>N,
    DIMRED = 0;
    fprintf(2,'flag\n');
    
end

%% ------------------------------------------------------------------------------

if ~isfield(options, 'w_true')
    GROUND_TRUTH = 0;
else
    GROUND_TRUTH = 1;
    w_true = options.w_true;
end



if M~=m,
    error('Sizes do not match')
end

y2sum = sum(y(:).^2)-mean(y).^2;
scale2 = y2sum / M;
% scale2 = y2sum / M/N;
% do not let it go very high, o/w the initial noise variance is estimated too low
scale2 = min([1e3,scale2]); 

scale = sqrt(scale2);


switch options.init,
    case 'ml'
        
        w = A \ y; 

        Sigma_w = zeros(N);
        z = ones(N,1);
        
        if strcmp(options.beta_init, 'auto')
            beta = 1./scale2;
        else
            beta = options.beta_init;
        end
    
    
    case 'rand'
        w = sqrt(scale)*randn(N,1);
        Sigma_w = zeros(N);
        z = scale*ones(N,1);
        
        beta = 1./scale2; 
        
end

% Precompute some expressions
AtA = A'*A;
Aty = A'*y;

%% Start iterations
tic

for it = 1:MAXITER,
    old_w = w;
    
    
    %% update X
    
    Z = diag(z);
    Zinv = diag(1./z);
    
    %% w step
    if cov_met == 1,
        Sigma_w = (beta*AtA + Z )^(-1);
        w = beta*Sigma_w*Aty;
            
    elseif cov_met == 2,
        Sigma_w = Zinv - Zinv*A'*(1/beta*eye(M)+A*Zinv*A')^(-1)*A*Zinv;
        w = beta*Sigma_w*Aty;
        
    elseif cov_met == 3,
        w = (beta*AtA + Z )\(beta*Aty);
%         w = Zinv*(eye(Neff) - A'*((1/beta*eye(M)+A*Zinv*A')\(A*Zinv)))*(beta*Aty);
        Sigma_w = diag(1./diag(beta*AtA + Z ));
    
    elseif cov_met == 4,
        w = (beta*AtA + Z )\(beta*Aty);
        Sigma_w = sparse(Neff,Neff);
    end
    
   
    
    %----------------------------------------------------------------------------------------
    %% estimate z
	
    dSigma_w = diag(Sigma_w);
    z = zeros(Neff,1);
    
    if INDIVIDUAL_GROUPS, % If all groups have size 1, for efficiency
        wg2 = w.^2 + dSigma_w;
        dg = 1;
        switch dist,
            case 'McKay'
                if strcmp(pt_estimate,'Mean')
                    z_g = (sqrt(a)./sqrt(wg2+eps)).*besselk(lambda-dg/2-1,sqrt(a).*sqrt(wg2)+eps)./besselk(lambda-dg/2,sqrt(a).*sqrt(wg2)+eps);
                    inv_z_g = (sqrt(wg2)./sqrt(a)).*besselk(lambda-dg/2+1,sqrt(a).*sqrt(wg2)+eps)./besselk(lambda-dg/2,sqrt(a).*sqrt(wg2)+eps);
                
                    if any(isnan(z_g))
                        idx = find(isnan(z_g));
                        z_g(idx) = sqrt(a(idx))./sqrt(wg2(idx)+eps);
                        inv_z_g(idx) = (sqrt(wg2(idx))./sqrt(a(idx)));
                    end
                         
                         
                elseif strcmp(pt_estimate,'Mode')
                    z_g = ( (dg/2+1-lambda)+sqrt((dg/2+1-lambda)^2 + a.*wg2) ) ./ wg2;
                    inv_z_g = 1./z_g;
                end
                
                
            case 'MLaplace'
                %lambda = (d+1)/2
                
                if strcmp(pt_estimate,'Mean')
                    z_g = sqrt(a)./sqrt(wg2+eps);
                    inv_z_g = sqrt(wg2)./sqrt(a).*besselk(3/2,sqrt(a).*sqrt(wg2)+eps)./besselk(1/2,sqrt(a).*sqrt(wg2)+eps);
                    if any(isnan(z_g))
                        idx = find(isnan(z_g));
                        z_g(idx) = sqrt(a(idx))./sqrt(wg2(idx)+eps);
                        inv_z_g(idx) = (sqrt(wg2(idx))./sqrt(a(idx)));
                    end
                    
                elseif strcmp(pt_estimate,'Mode')
                    z_g = ( (1/2)+sqrt(1/4 + a.*wg2) ) ./ wg2;
                    inv_z_g = 1./z_g;
                end
                
                
            case 'MStudent'
                if strcmp(pt_estimate,'Mean')
                    z_g = (dg/2-lambda)./(wg2/2 + b/2);
                elseif strcmp(pt_estimate,'Mode')
                    z_g = 2*(dg/2+1-lambda)./(b+wg2);
                end
                
            case 'Jeffreys'
                if strcmp(pt_estimate,'Mean')
                    z_g = dg./(wg2 + eps);
                elseif strcmp(pt_estimate,'Mode')
                    z_g = (dg+2)./(wg2+eps);
                end
        end
        z = z_g; 
        
    else
    
        for i=1:G
            wg2 = sum( w(grouping{i}).^2 + dSigma_w(grouping{i}) );
            dg = length(grouping{i});
            
            switch dist,
                case 'McKay'
                    if strcmp(pt_estimate,'Mean')
                        z_g(i) = (sqrt(a(i))./sqrt(wg2+eps)).*besselk(lambda-dg/2-1,sqrt(a(i)).*sqrt(wg2)+eps)./besselk(lambda-dg/2,sqrt(a(i)).*sqrt(wg2)+eps);
                        inv_z_g(i) = (sqrt(wg2)./sqrt(a(i))).*besselk(lambda-dg/2+1,sqrt(a(i)).*sqrt(wg2)+eps)./besselk(lambda-dg/2,sqrt(a(i)).*sqrt(wg2)+eps);
                        
                         if isnan(z_g(i))
                            z_g(i) =  (sqrt(a(i)))./(sqrt(wg2)+eps);
                            inv_z_g(i) = (sqrt(wg2)./sqrt(a(i)));
                         end
                         
                    elseif strcmp(pt_estimate,'Mode')
                        z_g(i) = ( (dg/2+1-lambda)+sqrt((dg/2+1-lambda)^2 + a(i)*wg2) ) ./ wg2;
                        inv_z_g(i) = 1./z_g(i);
                    end
                
                   
                
                case 'MLaplace'
                    if strcmp(pt_estimate,'Mean')
                        z_g(i) =  (sqrt(a(i)))./(sqrt(wg2)+eps);
                        inv_z_g(i) = (sqrt(wg2)./sqrt(a(i))).*besselk(3/2,sqrt(a(i)).*sqrt(wg2)+eps)./besselk(1/2,sqrt(a(i)).*sqrt(wg2)+eps);
                        
                        if isnan(inv_z_g(i))
                            inv_z_g(i) =  sqrt(wg2)./sqrt(a(i));
                        end
                         
                    elseif strcmp(pt_estimate,'Mode')
                        z_g(i) = ( (1/2)+sqrt(1/4 + a(i)*wg2) ) ./ (wg2+eps);
                        inv_z_g(i) = 1./z_g(i);
                    end
                    
                case 'MStudent'
                    if strcmp(pt_estimate,'Mean')
                        z_g(i) = 2*(dg/2-lambda)./(wg2 + b(i));
                    elseif strcmp(pt_estimate,'Mode')
                        z_g(i) = 2*(dg/2-lambda+1)./(wg2 + b(i));    
                    end
                    
                case 'Jeffreys'
                    if strcmp(pt_estimate,'Mean')
                        z_g(i) = dg./(wg2+eps);
                    elseif strcmp(pt_estimate,'Mode')
                        z_g(i) = (dg+2)./(wg2+eps);
                    end
            end
            
            
            % Implemented this way for overlapping groups
            z(grouping{i}) = z(grouping{i}) + z_g(i);
        
        end
    end
    
    z(isnan(z)) = 1e16;
    
    %% Estimate hyperparameters
    if ESTIMATE_HYPERPARAMS,
        switch dist,
            case 'McKay'
                if G == Neff
                    if strcmp(pt_estimate,'Mean')
                        a = (k_a+lambda)./(theta_a + inv_z_g/2);
                    elseif strcmp(pt_estimate,'Mode')
                        a = (k_a+lambda-1)./(theta_a + inv_z_g/2);
                    end
                else
                    for i=1:G,
                        if strcmp(pt_estimate,'Mean')
                            a(i) = (k_a+lambda)./(theta_a + inv_z_g(i)/2);
                        elseif strcmp(pt_estimate,'Mode')
                            a(i) = (k_a+lambda-1)./(theta_a + inv_z_g(i)/2);
                        end
                    end
                end
                
            case 'MLaplace'
                
                if G == Neff
                    if strcmp(pt_estimate,'Mean')
                        a = (k_a+1)./(theta_a + inv_z_g/2);
                    elseif strcmp(pt_estimate,'Mode')
                        a = (k_a)./(theta_a + inv_z_g/2);
                    end
                        
                else
                    for i=1:G,
                        dg = length(grouping{i});
                        if strcmp(pt_estimate,'Mean')
                            a(i) = (k_a+(dg+1)/2)./(theta_a + inv_z_g(i)/2);
                        elseif strcmp(pt_estimate,'Mode')
                            a(i) = (k_a+(dg+1)/2-1)./(theta_a + inv_z_g(i)/2);
                        end
                    end
                end
                
                
            case 'MStudent'
                
                if G == Neff
                    if strcmp(pt_estimate,'Mean')
                        b = (k_b - lambda)./ (theta_b + z_g/2);
                    elseif strcmp(pt_estimate,'Mode')
                        b = (k_b - lambda-1)./ (theta_b + z_g/2);
                    end
                else
                    for i=1:G,
                        if strcmp(pt_estimate,'Mean')
                            b(i) = (k_b - lambda)./ (theta_b + z_g(i)/2);
                        elseif strcmp(pt_estimate,'Mode')
                            b(i) = (k_b - lambda-1)./ (theta_a + z_g(i)/2);
                        end
                    end
                end
        end
    end
    
    
    
    
    %% update beta
    if  UPDATE_BETA & it>= UPDATE_BETA_init,
        err = sum(sum( abs(y - A*w).^2 ) ) + trace(AtA*Sigma_w) ;
        beta = (M + k_beta0)/(err+theta_beta0);
    end
    
   
    %% Prune irrelevant dimensions?
    if DIMRED,
        if INDIVIDUAL_GROUPS,
            
            
            if sum(find(z > DIMRED_THR)),
                indices = find(z <= DIMRED_THR);
                
                a = a(indices);
                b = b(indices);
                z = z(indices);
                w = w(indices);
                old_w = old_w(indices);
                A = A(:, indices);
                used_indices = used_indices(indices);
                
                AtA = A'*A;
                Aty = A'*y;
                [M G] = size(A);
                Neff = G;
                
            end
            
        else
            if sum(find(z_g > DIMRED_THR)),
                g_indices = find(z_g <= DIMRED_THR);
                
                
                indices = [];
                for i=1:length(g_indices)
                    N1 = length(indices);
                    N2 = length( grouping{g_indices(i)} );
                    indices = [indices; grouping{g_indices(i)}];
                    grouping{g_indices(i)} = [N1+1:N1+N2]';
                    
                end
                
                grouping = grouping(g_indices);
                G = length(g_indices);
                z_g = z_g(g_indices);
                
                a = a(g_indices);
                b = b(g_indices);
                z = z(indices);
                w = w(indices);
                old_w = old_w(indices);
                A = A(:, indices);
                used_indices = used_indices(indices);
                
                AtA = A'*A;
                Aty = A'*y;
                [M Neff] = size(A);
                
                
            end
        end
        
        
        
    end
    
    %%
    wconv = norm(old_w-w,'fro')/norm(old_w,'fro');
    
    if verbose,
        fprintf('max z = %g\n', max(z));
        
        if GROUND_TRUTH,
            werr = norm( w - w_true , 'fro') / norm( w_true , 'fro');
            fprintf('it %d: Error = %g, wconv = %g, beta = %g\n', it, werr, wconv,  beta);
        else
            fprintf('it %d: wconv = %g, beta = %g\n', it, wconv,  beta);
        end
    end
        
    if it > 50 & wconv < conv_thr 
        break;
    end
    
    if isnan(wconv),
        break;
    end
    
end
t = toc;


if DIMRED,
    temp = zeros(N,1);
    temp(used_indices) = w;
    w = temp;
    
    temp = zeros(N,1);
    temp(used_indices) = z;
    z = temp;
    
end

varargout{1} = it;
varargout{2} = z;
varargout{3} = beta;
varargout{4} = Sigma_w;
varargout{5} = a;
varargout{6} = b;
varargout{7} = used_indices;


function [options] = populate_vars(options)


if isempty(options)
    options.init                 = 'ml';
    options.verbose              = 0;
    options.lambda               = 1;
    options.a0                   = 1;
    options.b0                   = 1;
    options.k_a                  = 1e-6;
    options.theta_a              = 1e-6;
    options.k_b                  = 1e-6;
    options.theta_b              = 1e-6;
    options.k_beta0              = 1e-6;
    options.theta_beta0          = 1e-6;
    options.MAXITER              = 100;
    options.DIMRED               = 0;
    options.UPDATE_BETA          = 1;
    options.UPDATE_BETA_init     = 1;
    options.thr                  = 1e-10;
    options.dist                 = 'Jeffreys';
    options.pt_estimate          = 'Mean';
    options.cov_met              = 1;
    options.beta_init            ='auto';
    options.ESTIMATE_HYPERPARAMS = 1;
    options.conv_thr             = 1e-7;
    
else
    
    if ~isfield(options, 'init')
        options.init   = 'ml';
    end
    
    if ~isfield(options, 'verbose')
        options.verbose   = 0;
    end
    
    if ~isfield(options, 'lambda')
        options.lambda   = 1;
    end
    
    if ~isfield(options, 'a0')
        options.a0   = 1e-6;
    end
    
    if ~isfield(options, 'b0')
        options.b0   = 1e-6;
    end
    
    if ~isfield(options, 'k_a')
        options.k_a   = 1e-6;
    end
    
    if ~isfield(options, 'k_b')
        options.k_b   = 1e-6;
    end
    
    if ~isfield(options, 'theta_a')
        options.theta_a   = 1e-6;
    end
    
    if ~isfield(options, 'theta_b')
        options.theta_b   = 1e-6;
    end
    
    if ~isfield(options, 'k_beta0')
        options.k_beta0   = 1e-6;
    end
    
    if ~isfield(options, 'theta_beta0')
        options.theta_beta0   = 1e-6;
    end
    
    if ~isfield(options, 'MAXITER')
        options.MAXITER   = 100;
    end
    
    if ~isfield(options, 'DIMRED')
        options.DIMRED   = 0;
    end
    
    if ~isfield(options, 'UPDATE_BETA')
        options.UPDATE_BETA = 1;
    end
    
    if ~isfield(options, 'UPDATE_BETA_init')
        options.UPDATE_BETA_init = 1;
    end
    
    if ~isfield(options, 'ESTIMATE_HYPERPARAMS')
        options.ESTIMATE_HYPERPARAMS = 1;
    end
    
    if ~isfield(options,'dist')
        options.dist = 'Jeffreys';
    end
    
    if ~isfield(options,'cov_met')
        options.cov_met = 1;
    end
    
    if ~isfield(options,'pt_estimate')
        options.pt_estimate = 'Mean';
    end
    
    if ~isfield(options, 'beta_init')
        options.beta_init   = 'auto';
    end
    
   
    if ~isfield(options, 'conv_thr')
        options.conv_thr = 1e-7;
    end
    
    
    
end




