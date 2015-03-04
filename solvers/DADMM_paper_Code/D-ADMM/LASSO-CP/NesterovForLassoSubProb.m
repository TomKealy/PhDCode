function [xp, vars_prob] = NesterovForLassoSubProb(p, v, c, X, vars_prob)

% [xp, vars_prob] = NesterovForLassoSubProb(p, v, c, X, vars_prob)
%
% Solves the problem
%
%  minimize   phi(lambda) + (b/P + v)'*lambda + sig_P*t + c*||lambda||^2
%  subject to ||lambda|| <= t
%
% with variables x (vector of size n) and t (scalar). The function phi is
% given by
%
%    phi(lambda) = -inf_x ||x||_1 + (Ap'*lambda)'*x + delta/2*||x||^2.
%
% The gradient of the objective of the optimization problem is Lipschitz
% continuous and the Lipschitz constant is bounded by sigma_max(Ap)^2/delta
% + 2*c, where sigma_max(Ap) is the largest singular value of the matrix
% Ap. We use this information for solving this problem with Nesterov's
% method.


MAX_ITER = 5000;       % Maximum number of iterations of Nesterov's method
eps_stop = 1e-5;       % Stopping criterion parameter

% =========================================================================
% Input from vars_prob

n = vars_prob.n;       % Dimension of the primal variable (x in the paper)
max_sing_val_sq = vars_prob.max_sing_val_sq{p};  % Maximum singular value 
                                                 % squared
                                                 
one_over_P_b = vars_prob.one_over_P_b;           % (1/P)*b
sig_P = vars_prob.sig_P;                         % sigma/P
delta = vars_prob.delta;

previous_t = vars_prob.previous_t{p};

previous_lambda = X{p};
m =length(v);
% =========================================================================

% =========================================================================
% Pre-computations
Lipschitz_inv = 1/( max_sing_val_sq/delta + 2*c );

v_1_P_b = v + one_over_P_b;                      % v + (1/P)*b
zeros_vec = zeros(n,1);
% =========================================================================

% =========================================================================
% Initialization
eta = previous_lambda;
r = previous_t;

lambda = eta;
t = r;

k = 1;
prev_cost = Inf;
% =========================================================================

for iter = 1 : MAX_ITER
    
    factor = (k-1)/(k+2);
    
    % Compute the gradient of the cost function at (eta,r)
    [grad_phi, x_eta] = compute_grad_phi(eta, p, zeros_vec, vars_prob);
    grad_eta = [ grad_phi + v_1_P_b + 2*c*eta ; sig_P];
    
    new_point = [eta ; r] - Lipschitz_inv*grad_eta;
    
    % Project the new_point onto the second-order cone (or ice cream cone)
    lambda_prev = lambda;
    t_prev = t;
    [lambda, t] = lets_have_ice_cream(m, new_point);
            
    % ===============================================
    % Stopping criterion
    [grad_phi, x_lambda] = compute_grad_phi(lambda, p, zeros_vec, vars_prob);
    grad_lambda = [ grad_phi + v_1_P_b + 2*c*lambda ; sig_P];
    [proj_lambda, proj_t] = lets_have_ice_cream(m, [lambda;r]-grad_lambda);
    if max(proj_lambda - lambda) < eps_stop
        break;
    end
    % ===============================================
        
    eta = lambda + factor*(lambda - lambda_prev);
    r = t + factor*(t - t_prev);
    
%     % Adaptative restart
%     cost = -norm(x_lambda,1) + lambda'*grad_phi - (delta/2)*norm(x_lambda)^2 ...
%         + (v + one_over_P*b)'*lambda + sig_P*proj_t + c*norm(lambda)^2;
%     if cost > prev_cost
%         k = 0;
%     end
    
    k = k + 1;
    
end

if k == MAX_ITER
    fprintf('Warning: MAX_ITER reached in Nesterov''s method\n');
end

% =========================================================================
% Output

% For vars_prob
vars_prob.x_opt{p} = x_lambda;
vars_prob.previous_t{p} = t;

% For xp
xp = lambda;

% =========================================================================

end


function [grad, x_eta] = compute_grad_phi(eta, p, zero_vec, vars_prob)

% Computes the gradient of the function phi at the point eta. It returns
% the gradient and x_eta, the optimal value of x for this eta.

    Ap = vars_prob.A{p};
    delta = vars_prob.delta;
    
    x_eta = zero_vec;
    
    v_aux = Ap'*eta;
    
    ind_minus = v_aux < -1;
    ind_plus = v_aux > 1;
    
    x_eta(ind_minus) = -(v_aux(ind_minus) + 1)/delta;
    x_eta(ind_plus)  = -(v_aux(ind_plus) - 1)/delta;
    
    grad = Ap(:,ind_minus)*(-x_eta(ind_minus)) + Ap(:,ind_plus)*(-x_eta(ind_plus));
end


function [lambda, t] = lets_have_ice_cream(n, new_point)

% Projects new_point of dimensions n+1 onto the ice-cream cone.

x = new_point(1:n);
t = new_point(n+1);

norm_x = norm(x);

if norm_x <= t
    lambda = x;
    return;
end

if -norm_x < t && t < norm_x
    f = (t+norm_x)/2;
    lambda = f*x/norm_x;
    t = f;
    return;
end

% Else

lambda = zeros(n,1);
t = 0;

end







