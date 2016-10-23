clear all;

randn('seed', 0);
rand('seed',0);

n=2000;
m=500;
L=n;

positions = randi(L,[1,10]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 10;
Eb_N0_dB = 1;
S = randn(m,L);
A = S;
A = A/norm(A);
sigma = 10^(-Eb_N0_dB\20);
eta = randn(1,m);
noise_sum = sum(eta);
b = A*Tx_psd' + sigma*eta';

x0 = Tx_psd';

lambda = sqrt(sigma)*sqrt(2*log(n));

max_eig = max(abs(eig(A'*A)));

rho = nthroot(1/max_eig, 3);

max_iter = 500;

[x, history] = lasso(A, b, lambda, rho, 1.0);

[z, cost] = AMP(A, b, 2.0, x0, max_iter);

figure
semilogy(1:max_iter, cost, 'b', 1:max_iter, history.objval, 'r')
legend('AMP', 'ADMM');

figure
plot(1:n, z, 'b', 1:n, x, 'r')
legend('AMP', 'ADMM');