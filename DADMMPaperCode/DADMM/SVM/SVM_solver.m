function [x,vars_prob] = SVM_solver(p, v, c, X, vars_prob)

% [x,vars_prob] = SVM_solver(p, v, c, X, vars_prob)
%
% Solves the optimization problem
% 
%  minimize  1/(2*P)||s||^2 + beta*ones(1,m)*z + v'*[s;r] + c*||s||^2 
%                                                                   +c*r^2
%  subject to  y_k*(s'*x_k - r) >= 1 - z_k,  k = 1,...,m
%              z >= 0
%
% where x = (s,r) in R^n x R is the variable. s is the normal vector to an 
% hyperplane and r is a constant.
%
% The inputs p, v, c, and X are explained in the D-ADMM documentation. 
%
%

% *************************************************************************
% COMMENT THE NEXT CODE, IF YOU WANT THE BARZILAI-BORWEIN SOLVER
% *************************************************************************
% =========================================================================
% Initialization

% From struct vars_prob
A_qp = vars_prob.A_qp{p};
b_qp = vars_prob.b_qp;
Xip = vars_prob.Xi_p{p};
beta = vars_prob.beta;
P = vars_prob.P;
m = vars_prob.m;
n = vars_prob.n;
% =========================================================================


Hp_diag = [(1/P + 2*c)*ones(n,1) ; 2*c ; zeros(m,1)];
fp = [v ; beta*ones(m,1)];

x0 = [X{p} ; Xip];

opt_qp = optimset('Display', 'off', 'LargeScale', 'off');
x_qp = quadprog(diag(Hp_diag), fp, A_qp, b_qp, [], [], [], [], x0, opt_qp);

x = x_qp(1:n+1);

Xip = x_qp(n+2:end);
vars_prob.Xi_p{p} = Xip;
% *************************************************************************
% *************************************************************************




% *************************************************************************
% UNCOMMENT THE NEXT CODE, IF YOU WANT THE BARZILAI-BORWEIN SOLVER
% *************************************************************************
% 
% % Initialization
% 
% % From vars_prob
% A_p = vars_prob.A_p{p};
% m = vars_prob.m;
% n = vars_prob.n;
% P = vars_prob.P;
% beta = vars_prob.beta;
% warm_z = vars_prob.warm_z{p};
% 
% 
% % Pre-computed stuff
% Q = diag([(1/P + 2*c)*ones(n,1) ; 2*c]);
% Q_inv = diag([(P/(1+2*c*P))*ones(n,1) ; 1/(2*c)]);
% 
% AQA_inv = inv(A_p*Q_inv*A_p');
% 
% term_A = (Q_inv*A_p'*AQA_inv*A_p - eye(n+1))*Q_inv*v;
% term_B = Q_inv*A_p'*AQA_inv;
% term_C = term_A - term_B*ones(m,1);
% 
% aux1 = term_B'*Q*term_B;
% 
% U = [aux1 , -aux1 ; -aux1 , aux1];
% 
% QCv = Q*term_C+v;
% u1 = -term_B'*QCv;
% u2 = term_B'*QCv + beta*ones(m,1);
% 
% u = [u1 ; u2];
% 
% z = warm_z;
% Uz = U*z;
% grad_prev = Uz + u;
% f_prev = 0.5*z'*Uz + u'*z;
% 
% vars_BB = struct('f_val', {@SVM_BB_function}, ...
%     'grad', {@SVM_BB_gradient}, ...
%     'proj', {@SVM_BB_projection}, ...
%     'U', {U}, ...
%     'u', {u}, ...
%     'z_prev', {warm_z}, ...
%     'f_prev', {f_prev}, ...
%     'grad_prev', {grad_prev} ...
%     );
% 
% [z_opt, vars_BB] = SPG_BB(2*m, warm_z, vars_BB);
% 
% x = term_C + term_B*(z_opt(m+1:end) - z_opt(1:m));
% 
% vars_prob.warm_z{p} = z_opt;
% 
% *************************************************************************
% *************************************************************************