function [ w ] = islasso_SDPT3( XtX_over_n, XtY_over_n, G, D, IDX, lambda, params )

% ISLASSO SOLVED BY SECOND ORDER CONE PROGRAM - SDPT3
%
% --- we have tested this code with SDPT3 version 4.0(beta) (http://www.math.nus.edu.sg/~mattohkc/sdpt3.html) ---
%
% islasso_SDPT3( XtY_over_n, XtX_over_n, G, D, lambda, params )
%
% INPUTS:
%
% XtY_over_n  : X'*Y/n, correlations between design matrix and output vector Y (size p x 1)
% XtX_over_n  : X'*X /n, covariance matrix (size p x p)
% G           : group (logical) matrix (size p x ng, with ng the number of groups)
%               G(i,j) = true if variable i is in group j, false otherwise
% D           : weight matrix (same size as G)
%               D(i,j) = weight of the variable i in the group j if i is in j, 0 otherwise
% IDX         : index matrice that contains the sorting of G, per direction and size
%               (size nd x 2, with nd the number of directions in G; for instance 4 for the rectangular groups)
%                More details in the description of the function get_groups
% lambda      : regularization parameter
% params      : structure that contains additional paramaters (see the fields below)
%
% OUTPUTS:
% 
% w           : loading vector
%



%==========================================================================
if isfield(params,'threshold'), % Threshold for the variable selection step 
    threshold = params.threshold;
else
    threshold = 1e-8;
end
%==========================================================================

support = true(length(XtY_over_n),1);

for directions=1:size(IDX,1),

    tmp_idx = IDX(directions,1):IDX(directions,2);

    w  = slasso_SDPT3( XtX_over_n, XtY_over_n, G(:,tmp_idx), D(:,tmp_idx), lambda, params );

    support = support & ( abs(w) > threshold );

end

w(~support) = 0;
w(support)  = ( XtX_over_n(support,support) )\( XtY_over_n(support) );
