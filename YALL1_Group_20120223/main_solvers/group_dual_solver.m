function [x,Out] =  group_dual_solver_20110627(A,At,invAAt,b,groups,g,opts)
%
%  Dual-based ADM solver for the non-overlapping group-sparse problem:
%
%  minimize     sum_i w_i*||x_{gi}||_2
%  subject to   Ax = b,
%               x >= 0 (optional),
%
%  where gi is the index set of the i-th group, w_i is the 
%  weight for the i-th group, and the groups do not overlap.
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
m = size(b,1);
n = length(groups);
if isfield(opts,'xInit')
    x = opts.xInit;
else
    x = zeros(n,1);
end
y = zeros(m,1);
Aty = zeros(n,1);

%--------------------------------------------------------------
% Main iterations
%--------------------------------------------------------------
for iter = 1:maxit
    
    % Solve z-subproblem by group-wise projection
    z = proj_group(Aty+x/beta,weights,nonneg,groups,g);
    
    % Solve y-subproblem (convex quadratic)
    [y,Aty] = linsolve_dual(A,At,b,x,z,beta,opts,invAAt,y,Aty);
    
    % Check the stopping criterion
    r = z-Aty;
    xd = gamma*beta*r;    
    if iter > 20 && ~mod(iter,5)    % check every 5 iterations starting at 25
        crit_met = norm(xd) < tol*norm(x);
        
        % --- put your stopping condition here ---
        % crit_met = ......

        if crit_met
            Out.status = 'Stopping criterion met';
            break; 
        end
    end
    % Update x 
    x = x - xd;
    
    % Continuation on the penalty parameter
    if opts.continuation
        if iter > 1
            beta = continuation_dual(beta,r,rp,opts);
        end
        rp = r;
    end
    
end

% Outputs
Out.iter = iter;
if iter==maxit, Out.status = 'Maximum # of iterations reached'; end

end

