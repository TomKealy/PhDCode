%
% Code for generating time-frequency separation example
%
clear all;

n = 2^10;

k = 5;

I = eye(n);
F = RSTMatrix(n); % orthonormal

A = [I F];

x = zeros(2*n,1);
p = randperm(n);
x(p(1:k)) = rand(k,1);
p = randperm(n);
x(n+p(1:k)) = rand(k,1);

y = A*x;

tic
xhat_ist = SolveIST(A, y);
toc

yhat_ist = A*xhat_ist; % reconstruction from components solution
figure; plot(y,'-.r','LineWidth',2); hold on; plot(yhat_ist); axis([0 n -1/4 1]); print('-depsc', 'bp-ist-signal.eps');
figure; plot(x(1:n),'-.r','LineWidth',2); hold on; plot(xhat_ist(1:n)); axis([0 n 0 1]); print('-depsc', 'bp-ist-dirac.eps');
figure; plot(x(n+1:2*n),'-.r', 'LineWidth',2); hold on; plot(xhat_ist(n+1:2*n)); axis([0 n 0 1]); print('-depsc', 'bp-ist-fourier.eps');

rmse_error_ist = sqrt(norm(x - xhat_ist)/length(x)) % error
rmse_resid_ist = sqrt(norm(y - yhat_ist)/length(y)) % residual

% compare time and error with fast version of LARS using rank one updates
tic
[xhat_lars, numIters, activationHist] = FastLars(A, y);
toc

yhat_lars = A*xhat_lars; % reconstruction from components solution
figure; plot(y,'-.r','LineWidth',2); hold on; plot(yhat_lars); axis([0 n -1/4 1]); print('-depsc', 'bp-lars-signal.eps');
figure; plot(x(1:n),'-.r','LineWidth',2); hold on; plot(xhat_lars(1:n)); axis([0 n 0 1]); print('-depsc', 'bp-lars-dirac.eps');
figure; plot(x(n+1:2*n),'-.r', 'LineWidth',2); hold on; plot(xhat_lars(n+1:2*n)); axis([0 n 0 1]); print('-depsc', 'bp-lars-fourier.eps');
numIters
rmse_error_lars = sqrt(norm(x - xhat_lars)/length(x)) % error
rmse_resid_lars = sqrt(norm(y - yhat_lars)/length(y)) % residual

%
% Copyright (c) 2006. Iddo Drori
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
