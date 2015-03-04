function [xp, vars_prob] = BPDN_RP_Solver(p, v, c, X, vars_prob)

% [xp, vars_prob] = solver_BPDN(p, v, c, X, vars_prob)
%
% Solves the problem
%
%         minimize  fp(x) + v'*x + c*||x||^2 ,                 (1)
%            x
%
% where fp(x) = 0.5*||Ap*x - bp||^2 + (beta/P)*||x||_1 is the individual's
% node function.
%
% Doing algebraic manipulations, it can be shown that (1) is equivalent to
%
%         minimize  0.5*||B*x - r||^2 + (beta/P)*||x||_1       (2)
%            x
%
% where B = [Ap ; sqrt(2*c)*eye(n)], and r = [bp ; -1/sqrt(2*c)*v]. We 
% solve (2) with the GPRS_BB algorithm. We will use the previous solution 
% as a warm-start.

A_full = vars_prob.A_BPDN;
b_full = vars_prob.b_BPDN;
m_p    = vars_prob.m_p; 
P      = vars_prob.P;
beta   = vars_prob.beta;
handler_GPSR = vars_prob.handler_GPSR;

Ap = A_full(1+(p-1)*m_p : p*m_p , :);
bp = b_full(1+(p-1)*m_p : p*m_p);

n = length(X{p});

%B = [Ap ; sqrt(2*c)*eye(n)];
r = [bp ; -(1/sqrt(2*c))*v];

B_handler  = @(x) [Ap*x ; sqrt(2*c)*x];
BT_handler = @(x) Ap'*x(1:m_p) + sqrt(2*c)*x(m_p+1:end);

xp = handler_GPSR(r, B_handler, beta/P, 'Initialization', X{p}, 'Verbose', 0, ...
    'StopCriterion', 3, 'ToleranceA', 1e-6, 'AT', BT_handler);






