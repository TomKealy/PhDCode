function [s1,x]=TCM_state_eq1(s,u,Constellation) 
% State equation for the TCM encoder in Fig. P9.9.1(a)
% Input:  s= State, u= Input
% Output: s1= Next state, x= Output
s1 = [u(2) s(1)];
x= Constellation(bin2deci([u(1) rem(u(2)+s(2),2) s(1)])+1);