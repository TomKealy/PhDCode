close all; clear all;

m = 500;
M = 1000;

data = csvread('tvws_data1.csv');

data = data(1:8000, 2);

plot(data)

A = randn(m,M);
F = LehmerMatrix(M);
[L U] = lu(F);
I = eye(M);
D = inv(L);

data = reshape(data, 8, 1000);

Y = A*data';

rho = 1/(max(abs(eig(A'*A)))); 
rho = nthroot(rho, 3);
lambda = sqrt(2*log(1000))/8;

[B, cost] = lasso_admm_lucy(Y, A*L, lambda, 0.5);

C = reshape(B, 1, 8000);

figure
plot(C)
title('B')