function [beta,invIpAAt] = continuation_primal_20110620(beta,r1,r1p,r2,r2p,invIpAAt,opts)
% 
%  Continuation on penalty parameters for primal solvers.
%  
%  Inputs: 
%         beta - penalty parameters
%           r1 - constraint violation z-x at the current iteration
%          r1p - constraint violation z-x at the previous iteration
%           r2 - constraint violation Ax-b at the current iteration
%          r2p - constraint violation Ax-b at the previous iteration
%     invIpAAt - a function handle for (beta1*eye(m)+beta2*(A*A'))\x
%         opts - a structure with continuation parameters
%
%  Outputs:
%         beta - updated penalty parameters
%     invIpAAt - updated function handle
%
%--------------------------------------------------------------------------

% Get continuation parameters
contpar    = opts.contpar;
contfactor = opts.contfactor;

% Perform continuation
if opts.nonorth && opts.exact
    % In this case, penalty parameters beta(1) and beta(2) are updated
    % simultaneously for simple update of invIpAAt
    if norm(r1(:))>contpar*norm(r1p(:)) && norm(r2(:))>contpar*norm(r2p(:))
        beta = contfactor*beta;
        invIpAAt = @(x) invIpAAt(x)/contfactor;
    end
else
    if norm(r1(:)) > contpar*norm(r1p(:))
        beta(1) = contfactor*beta(1);
    end
    if norm(r2(:))> contpar*norm(r2p(:))
        beta(2) = contfactor*beta(2);
    end
end

end