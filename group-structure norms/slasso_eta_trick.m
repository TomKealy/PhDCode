function [ w ] = slasso_eta_trick( X, Y, XtY_over_n, XtX_over_n, G, D, lambda, params )

% SLASSO SOLVED BY FIRST ORDER TECHNIQUE (see Appendix G in the paper)
%
% slasso_eta_trick( X, Y, XtY_over_n, XtX_over_n, G, D, lambda, params )
%
% INPUTS:
%
% X           : design matrix (size n x p, i.e., n row data of dimension p)
% Y           : output vector (size n x 1)
% XtY_over_n  : X'*Y/n, correlations between design matrix and output vector Y (size p x 1)
% XtX_over_n  : X'*X /n, covariance matrix (size p x p)
% G           : group (logical) matrix (size p x ng, with ng the number of groups)
%               G(i,j) = true if variable i is in group j, false otherwise
% D           : weight matrix (same size as G)
%               D(i,j) = weight of the variable i in the group j if i is in j, 0 otherwise
% lambda      : regularization parameter
% params      : structure that contains additional paramaters (see the fields below)
%
% OUTPUTS:
% 
% w           : loading vector
%
%   Copyright (c) 2009 Rodolphe Jenatton. All rights reserved.

YtY_over_n = Y'*Y / length(Y) ;

%==========================================================================
% Process the params
%==========================================================================
if isfield(params,'epsilon'), % Smoothing parameter for the closed-form update (see Appendix G)
    epsilon = params.epsilon;
else
    epsilon        = 1e-6;
    params.epsilon = epsilon;
end
%--------------------------------------------------------------------------
if isfield(params,'max_it'),% Maximum number of iterations
    max_it = params.max_it;
else
    max_it        = 300; 
    params.max_it = max_it;
end
%--------------------------------------------------------------------------
if isfield(params,'min_delta_cost'), % Minimum relative descrease between two iterations
    min_delta_cost = params.min_delta_cost;
else
    min_delta_cost        = 1e-6;
    params.min_delta_cost = min_delta_cost;
end
%--------------------------------------------------------------------------
if isfield(params,'w'),% Loading vector for a potentiel warm-restart
    w = params.w;
else
    w = zeros(size(X,2),1);
end
%--------------------------------------------------------------------------

t          = 1;
t0         = 10;
delta_cost = Inf;
 
while ( t <= max_it) && ( delta_cost > min_delta_cost), 
    
    %======================================================================
    % Update of Zeta
    
    Zeta = mexEta_update_for_regression( w, G, D, epsilon/t );
    
    %======================================================================
    % Update of w
    
    w    = mexW_Regularized_update( X, Y, Zeta, lambda );
        
    %======================================================================
    % Compute cost function
    
    if ( mod(t,t0) == 0 ),
        
      cost = YtY_over_n / 2 + ( w'*XtX_over_n*w )/2 - XtY_over_n'*w  + lambda * get_norm_evaluation( w, G, D );
      
      if ( t > t0 ),
          
        delta_cost = 100*(past_cost - cost)/past_cost;
        
        if (delta_cost < 0) && ( delta_cost > -1e-6 ), delta_cost = -delta_cost; end;
        
      end
      
      past_cost  = cost;
    
      fprintf( 'Iteration %4.3d, cost function : %8.6f, cost variation: %8.6f \n',t, cost, delta_cost); 
    
    end
    %======================================================================
    t = t + 1;
    
end


end
