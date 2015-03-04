clear;
close all;

n = 512;    % size of signal
S = 28;     % Sparsity

support = randsample(n,S);
x0 = zeros(n,1); x0(support) = randn(S,1);
figure(1); hold on
plot(1:n,x0,'b*');

% Sampling
m = 64;
k=0:(n-1); t=0:(n-1);
F = exp(-j*2*pi*(k'*t)/n)/sqrt(n);  % Fourier matrix
freq = randsample(n,m);
A = [real(F(freq,:)); imag(F(freq,:))];

% Noisy projection
noise_bound = (1e-1)*rand();
noise = rand(size(A,1),1); noise = noise/norm(noise);
y = A*x0 + noise_bound*noise;

% l-1 minimization via l-1 Magic
addpath l1magic/Optimization
x1=A\y;     % Initial value via l-2 minimization
figure(2);
plot(1:n, x1,'b*');

x = l1qc_logbarrier(x1,A,[],y,noise_bound, 1e-3); % second-order cone programming
figure(1);
plot(1:n,x,'ro');