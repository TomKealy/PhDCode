function [x,vars_prob] = AverageConsensusSubProb(p, v, c, X, vars_prob)

% [x,vars_prob] = AverageConsensusSubProb(p, v, c, X, vars_prob)
%
% Solves the quadratic problem
% 
%             minimize    (x-theta_p)^2 + vx + cx^2      (1)  
%                x
%
% This function is used to solve the average consensus problem in a
% distributed way.
%
% The inputs p, v, c, and X are explained in the D-ADMM documentation. 

theta_p = vars_prob.theta(p);

x = (2*theta_p - v)/(2*(1+c));


