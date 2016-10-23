function [ w, primal_,dual_,slack_, info] = slasso_SDPT3( XtX_over_n, XtY_over_n, G, D, lambda, params )

% SLASSO SOLVED BY SECOND ORDER CONE PROGRAM - SDPT3
%
% --- we have tested this code with SDPT3 version 4.0(beta) (http://www.math.nus.edu.sg/~mattohkc/sdpt3.html) ---
%
% slasso_SDPT3( XtY_over_n, XtX_over_n, G, D, lambda, params )
%
% INPUTS:
%
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
% The other outputs gather the results of the optimization performed by
% SDPT3 (see documentation of SDPT3). This enables warm-restart.


%==========================================================================
% Params processing

if isfield(params,'constrained_mode'), % Solved SLASSO in a constrained formulation 
    constrained_mode = params.constrained_mode;
else
    constrained_mode = false;
end
%--------------------------------------------------------------------------
if isfield(params,'squared_penalization_mode'), % Solved SLASSO with the squared-penalization formulation 
    squared_penalization_mode = params.squared_penalization_mode;
else
    squared_penalization_mode = false;
end
%--------------------------------------------------------------------------
if isfield(params,'warm_start') && ... % Variables required by SDPT3 solver to perform warm-restart
   params.warm_start            && ...
   isfield(params,'primal')     && ...
   isfield(params,'dual')       && ...
   isfield(params,'slack'),
    
    warm_start = true;
    
    primal = params.primal;
    dual   = params.dual;
    slack  = params.slack;
   
else
    warm_start = false;
end
%--------------------------------------------------------------------------
if isfield(params,'lambda_normalization'), % It might useful to normalize the regularization parameter
    lambda_normalization = params.lambda_normalization;
else
    lambda_normalization = 1;
end
%==========================================================================

p = length(XtY_over_n);

if constrained_mode,
    
    [ At, b, c, blk ] = get_sdpt3_chol_constrained_setting(XtX_over_n, XtY_over_n, G, D, lambda/lambda_normalization);
    
elseif squared_penalization_mode,
    
    [ At, b, c, blk ] = get_sdpt3_chol_sq_regularized_setting(XtX_over_n, XtY_over_n, G, D, lambda/lambda_normalization);
    
else
    
    [ At, b, c, blk ] = get_sdpt3_chol_regularized_setting(XtX_over_n, XtY_over_n, G, D, lambda/lambda_normalization);
    
end

OPTIONS.printlevel = 0;

if warm_start,
    [obj, primal_,dual_,slack_, info] = sqlp(blk, At, c, b, OPTIONS, primal, dual, slack);
else
    [obj, primal_,dual_,slack_, info] = sqlp(blk, At, c, b, OPTIONS);
end

w = dual_(end-p+1:end);
