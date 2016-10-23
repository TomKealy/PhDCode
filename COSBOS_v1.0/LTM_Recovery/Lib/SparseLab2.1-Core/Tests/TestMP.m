% Test the routine SolveMP
clear all;

n = 200;
d = 100;
k = 20;

A = MatrixEnsemble(d,n);
x = 100*[rand(k,1)-0.5*ones(k,1); zeros(n-k,1)];
y = A*x;

[coefMP2, i2, ah2] = SolveMP(A, y, n, 10000, 0, 0, 1);
    
%
% Copyright (c) 2006. Yaakov Tsaig
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
