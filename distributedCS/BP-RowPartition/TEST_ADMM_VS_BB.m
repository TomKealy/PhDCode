
m = 200;
n = 1000;

A = randn(m,n);
b = randn(m,1);
dual_var = {zeros(m,1)};

vars_prob = struct('A_BP', {A}, ...
    'b_BP', {b}, ...
    'dual_var', {dual_var} ...
    );

v = randn(n,1);
c = 10;

t_BB_aux = cputime;
[x_BB,vars_prob] = minimize_quad_prog_plus_l1_BB(1, v, c, 1, vars_prob);
t_BB = cputime - t_BB_aux;

% =========================================================================
% ADMM

P = 1;

A_pinv = pinv(A);

% *****************
% Parameters
MAX_ITER = 500;
rho = 1;
tau_rho = 10;
mu_rho = 2;

eps_prim = 1e-6;
eps_dual = 1e-5;
% *****************

lambda = zeros(n,1);

x_ADMM = zeros(n,1);
y_ADMM = zeros(n,1);

t_ADMM_aux = cputime;

for i = 1 : MAX_ITER
   
   % **********************************************************
   % x-minimization
   v_new = v + lambda - rho*y_ADMM;
   c_new = c + 0.5*rho;
   
   x_ADMM = zeros(n,1);
   
   ind_pos = (v_new < -(1/P));
   ind_neg = (v_new >  (1/P));
   
   x_ADMM(ind_pos) = -(1/(2*c_new))*( v_new(ind_pos) + (1/P) );
   x_ADMM(ind_neg) = -(1/(2*c_new))*( v_new(ind_neg) - (1/P) );
   % **********************************************************
   
   % **********************************************************
   % y-minimization
   
   y_prev = y_ADMM;
   
   point = (1/rho)*lambda + x_ADMM;   
   y_ADMM = point + A_pinv*(b - A*point);
   % **********************************************************
   
   % **********************************************************
   % lambda update
   
   r_prim = x_ADMM - y_ADMM;    % primal residual   
   lambda = lambda + rho*r_prim;
   
   s_dual = rho*(y_ADMM - y_prev);  % dual residual      
   % **********************************************************
   
   % **********************************************************
   % rho adjustment
   
   r_prim_norm = norm(r_prim);
   s_dual_norm = norm(s_dual);
   
   if r_prim_norm > tau_rho*s_dual_norm
       rho = mu_rho*rho;
   elseif s_dual_norm > tau_rho*r_prim_norm
       rho = rho/mu_rho;
   end
   % **********************************************************
   
   if r_prim_norm < eps_prim && s_dual_norm < eps_dual
       break;
   end    
end
t_ADMM = cputime - t_ADMM_aux;
% =========================================================================


cvx_begin
  variable x_cvx(n,1);
  minimize((1/P)*norm(x_cvx,1) + v'*x_cvx + c*square_pos(norm(x_cvx)));
  subject to
    A*x_cvx == b;
cvx_end

fprintf('Error BB   = %f\n', norm(x_BB - x_cvx)/norm(x_cvx));
fprintf('Error ADMM = %f\n', norm(x_ADMM - x_cvx)/norm(x_cvx));
fprintf('Time BB    = %f\n', t_BB);
fprintf('Time ADMM  = %f\n', t_ADMM);


