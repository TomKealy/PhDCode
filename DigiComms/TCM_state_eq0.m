function [s1,x]=TCM_state_eq0(s,u,Constellation) 
% State equation for the TCM_encoder in Fig. 9.14(a)
% Input:  s= State, u= Input, Constellation
% Output: s1= Next state, x= Output
s1= [s(3) rem([s(1)+u(1) s(2)+u(2)],2)];
x= Constellation(bin2deci([u(1) u(2) s(3)])+1);