function [x,Out] = overlap_primal_solver_20120219(A,At,b,groups,gmat,G,opts)
%
%  Primal-based ADM solver for the overlapping group-sparse problem:
%
%  minimize     sum_i w_i*||x_{gi}||_2
%  subject to   Ax = b,
%               x >= 0 (optional),
%
%  where gi is the index set of the i-th group, w_i is the weight 
%  for the i-th group, and the groups may overlap.
%
% ---------------------------------------------------------------------

% Get optional parameters
tol          = opts.tol;
weights      = opts.weights;
nonneg       = opts.nonneg;
maxit        = opts.maxit;
beta         = opts.beta;  % must be a 2-vector
gamma        = opts.gamma; % must be a 2-vector

% Initialization
m = size(b,1);
[N,n] = size(G);
if isfield(opts,'xInit')
    x = opts.xInit;
    Ax = A(x); Gx = G*x;
else
    x = zeros(n,1); Ax = zeros(m,1); Gx = zeros(N,1);
end
lambda1 = zeros(N,1);
lambda2 = zeros(m,1);

%--------------------------------------------------------------
% Main iterations
%--------------------------------------------------------------
for iter = 1:maxit
    
    % Solve z-subproblem by group-wise shrinkage
    z = shrink_group(Gx+lambda1/beta(1),weights/beta(1),nonneg,groups,gmat);
    
    % Solve x-subproblem (convex quadratic)
    xp = x;  % record the previous iterate
    [x,Ax,Gx] = linsolve_overlap_primal(A,At,G,b,z,beta,lambda1,lambda2,...
                                        xp,Ax,Gx);
    
    % Check the stopping criterion
    if iter > 20 && ~mod(iter,5)    % check every 5 iterations starting at 25
        crit_met = norm(x-xp) < tol*norm(xp);

        % --- put your stopping condition here ---
        % crit_met = ......

        if crit_met
            Out.status = 'Stopping criterion met';
            break; 
        end
    end
    
    % Update the multipliers
    r1 = z-Gx; r2 = Ax-b;
    lambda1 = lambda1-gamma(1)*beta(1)*r1;
    lambda2 = lambda2-gamma(2)*beta(2)*r2;
    
    % Continuation on the penalty parameters
    if opts.continuation
        if iter > 1
            [beta,~] = continuation_primal(beta,r1,r1p,r2,r2p,[],opts);
        end
        r1p = r1; r2p = r2;
    end
end

% Outputs
Out.iter = iter;
if iter==maxit, Out.status = 'Maximum # of iterations reached'; end

end


