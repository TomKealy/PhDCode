function beta = continuation_dual_20110620(beta,r,rp,opts)
% 
%  Continuation on the penalty parameter for dual solvers.
%  
%  Inputs: 
%         beta - penalty parameter
%            r - constraint violation z-A'y at the current iteration
%           rp - constraint violation z-A'y at the previous iteration
%         opts - a structure with continuation parameters
%
%  Outputs:
%         beta - updated penalty parameter
%
%--------------------------------------------------------------------------

% Get continuation parameters
contpar    = opts.contpar;
contfactor = opts.contfactor;

% Perform continuation
if norm(r(:)) > contpar*norm(rp(:))
    beta = contfactor*beta;
end

end