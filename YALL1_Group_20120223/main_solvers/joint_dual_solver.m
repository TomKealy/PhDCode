function [X,Out] =  joint_dual_solver_20110627(A,At,invAAt,B,opts)
%
%  Dual-based ADM solver for the jointly-sparse basis pursuit model:
%
%  minimize     sum_i w_i*||X(i,:)||_2
%  subject to   AX = B,
%               X >=0 (optional),
%
%  where X(i,:) is the i-th row of matrix X, and w_i is the
%  weight for the i-th row.
%
% -------------------------------------------------------------

% Get optional parameters
tol          = opts.tol;
weights      = opts.weights;
nonneg       = opts.nonneg;
maxit        = opts.maxit;
beta         = opts.beta;   % must be a scalar
gamma        = opts.gamma;  % must be a scalar

% Initialization
[m,L] = size(B);
n = length(At(zeros(m,1)));
if isfield(opts,'xInit')
    X = opts.xInit;
else
    X = zeros(n,L);
end
Y = zeros(m,L);
AtY = zeros(n,L);

%--------------------------------------------------------------
% Main iterations
%--------------------------------------------------------------
for iter = 1:maxit
    
    % Solve Z-subproblem by row-wise projection
    Z = proj_group(AtY+X/beta,weights,nonneg);
    
    % Solve Y-subproblem (convex quadratic)
    [Y,AtY] = linsolve_dual(A,At,B,X,Z,beta,opts,invAAt,Y,AtY);
    
    % Check the stopping criterion
    R = Z-AtY;
    Xd = gamma*beta*R;
    if iter > 20 && ~mod(iter,5)    % check every 5 iterations starting at 25
        crit_met = norm(Xd(:)) < tol*norm(X(:));

        % --- put your stopping condition here ---
        % crit_met = ......

        if crit_met;
            Out.status = 'Stopping criterion met';
            break;
        end
    end
    % Update X
    X = X - Xd;
    
    % Continuation on the penalty parameter
    if opts.continuation
        if iter > 1
            beta = continuation_dual(beta,R,Rp,opts);
        end
        Rp = R;
    end
    
% Outputs
Out.iter = iter;
if iter==maxit, Out.status = 'Maximum # of iterations reached'; end

end

