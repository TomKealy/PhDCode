% This program simulates a fast compressed sensing scheme.
% It generates an implicit random sampling matrix of size nxN. 
% It then approximately solves the CS problem
%   Min || alpha ||_1 s.t. Phi*alpha = Phi*alpha0
% using various solvers.
% The CS matrix is constructed as P*H*Q, where P,Q are random permutation
% matrices, and H is a fast Hadamard/Fourier operator.

n = 500;
N = 1000;
k = 30;

x0 = [randn(k,1); zeros(N-k,1)];

% Calculate y using the CS matrix
y = FastCSOperator(1,n,N,x0,1:N,N);

% Initialize threshold parameters
delta = n/N;
rho = k/n;
S = 10;
alpha_0 = delta*(1-rho)/S;
q = min((n-k)/k,0.5);

% Solve the CS problem 
[xhat, numIters] = SolveStOMP('FastCSOperator', y, N, 'FAR', alpha_0, S, 1);
[xhat2, numIters] = SolveLasso('FastCSOperator', y, N, 'lars');
[xhat3, numIters] = SolveOMP('FastCSOperator', y, N);
xhat4 = SolveBP('FastCSOperator', y, N);

figure;
subplot(5,1,1); plot(x0(1:100)); title('Original');
subplot(5,1,2); plot(xhat(1:100)); title(['ITSP, Error = ' num2str(norm(x0-xhat),3)]);
subplot(5,1,3); plot(xhat2(1:100)); title(['Lars, Error = ' num2str(norm(x0-xhat2),3)]);
subplot(5,1,4); plot(xhat3(1:100)); title(['OMP, Error = ' num2str(norm(x0-xhat3),3)]);
subplot(5,1,5); plot(xhat4(1:100)); title(['BP, Error = ' num2str(norm(x0-xhat4),3)]);
    
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
