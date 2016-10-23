function [X,Out] = joint_primal_solver(A,At,invIpAAt,B,opts)
%
%  Primal-based ADM solver for the jointly-sparse problem:
%
%  minimize     sum_i w_i*||X(i,:)||_2
%  subject to   AX = B,
%               X >= 0 (optional),
%
%  where X(i,:) is the i-th row of matrix X, and w_i is the
%  weight for the i-th row.
%
% -------------------------------------------------------------

% Get optional parameters
tol          = opts.tol;
weights      = opts.weights;
nonneg       = opts.nonneg;
nonorth      = opts.nonorth;
exact        = opts.exact;
maxit        = opts.maxit;
beta         = opts.beta;  % must be a 2-vector
gamma        = opts.gamma; % must be a 2-vector

% Initialization
[m,L] = size(B);
n = length(At(zeros(m,1)));
if isfield(opts,'xInit')
    X = opts.xInit;
    if nonorth && ~exact, AX = A(X); end
else
    X = zeros(n,L); AX = zeros(m,L);
end
Lambda1 = zeros(n,L);
Lambda2 = zeros(m,L);

%--------------------------------------------------------------
% Main iterations
%--------------------------------------------------------------
for iter = 1:maxit
    
    % Solve Z-subproblem by row-wise shrinkage
    Z = shrink_group(X+Lambda1/beta(1),weights/beta(1),nonneg);
    
    % Solve X-subproblem (convex quadratic)
    Xp = X;  % record the previous iterate
    [X,AX] = linsolve_primal(A,At,B,Z,beta,Lambda1,Lambda2,opts,invIpAAt,Xp,AX);
        
    % Check the stopping criterion
    if iter > 20 && ~mod(iter,5)    % check every 5 iterations starting at 25
        crit_met = norm(X(:)-Xp(:)) < tol*norm(Xp(:));

        % --- put your stopping condition here ---
        % crit_met = ......

        if crit_met
            Out.status = 'Stopping criterion met';
            break; 
        end
    end
    
    % Update the multipliers
    r1 = Z-X; r2 = AX-B;
    Lambda1 = Lambda1-gamma(1)*beta(1)*r1;
    Lambda2 = Lambda2-gamma(2)*beta(2)*r2;
    
    % Continuation on the penalty parameters
    if opts.continuation
        if iter > 1
            [beta,invIpAAt] = continuation_primal(beta,r1,r1p,r2,r2p,invIpAAt,opts);
        end
        r1p = r1; r2p = r2;
    end
end

% Outputs
Out.iter = iter;
if iter==maxit, Out.status = 'Maximum # of iterations reached'; end

end