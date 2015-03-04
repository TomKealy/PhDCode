n=200;
m=50;
L=n;

positions = randi(L,[1,5]);%generate random spikes for signal

Tx_psd = zeros(1,L); %Tx PSD
Tx_psd(positions) = 1;
Eb_N0_dB = 10;
S = randn(m,L);
A_BP = S;
sigma = 10^(-Eb_N0_dB\20);
eta = randn(1,m);
noise_sum = sum(eta);
b = A_BP*Tx_psd' + sigma*eta';

% For the groundtruth, we use the spgl1 solver

lambda = 0.01*norm(A_BP'*b, 'inf');

%solution = spgl1(A_BP, b, 0, lambda, []);

s = lasso(A_BP, b, lambda, 1, 1); 
