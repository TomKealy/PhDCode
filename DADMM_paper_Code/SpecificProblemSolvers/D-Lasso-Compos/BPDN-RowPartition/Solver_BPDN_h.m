function [xp, vars_prob] = Solver_BPDN_h(p, v, c, X, vars_prob)

% [xp, vars_prob] = Solver_BPDN_h(p, v, c, X, vars_prob)
%
% Solves the problem
%
%         minimize  hp(x) + v'*x + c*||x||^2 ,                 (1)
%            x
%
% where hp(x) = 0.5*||Ap*x - bp||^2.
%


Apbp = vars_prob.Apbp{p};
M_inv = vars_prob.M_inv{p};

vec = Apbp - v;
xp = M_inv*vec;






