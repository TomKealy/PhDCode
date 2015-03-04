function [xp, vars_prob] = Solver_BPDN_g(p, v, c, X, vars_prob)

% [xp, vars_prob] = Solver_BPDN_g(p, v, c, X, vars_prob)
%
% Solves the problem
%
%         minimize  gp(x) + v'*x + c*||x||^2 ,                 (1)
%            x
%
% where gp(x) = (beta/P)*||x||_1.

P      = vars_prob.P;
beta   = vars_prob.beta;

n = length(v);

v_t = (P/beta)*v;
c_t = (P/beta)*c;

xp = zeros(n,1);

ind_minus_1 = (v_t < -1);
ind_plus_1  = (v_t >  1);

xp(ind_minus_1) = -(v_t(ind_minus_1) + 1)/(2*c_t);
xp(ind_plus_1) = -(v_t(ind_plus_1) - 1)/(2*c_t);


