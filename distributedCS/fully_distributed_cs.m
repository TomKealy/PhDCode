clear all
close all

% W=6000 bin singal divided into L=500 subbands of bandwidth B=12 bins
% (n.b. all in Fourier domain).

W=600; %bandwidth of spectrum
L=201; %number of sub-bands
L0 = (L-1)/2; %re-labelling via shift
J = 50; %number of nodes
B=W/L; %bandwidth per subband
K = 0.02*L; % sparsity

positions = randi(L,[1,K]);%generate random spikes for signal

s = zeros(1,L);
s(positions) = 1;
figure
plot(s)

% fs = fft(s);
% figure
% plot(abs(fs))
% % figure
% % plot(imag(fs))

%Secondary User system

h = sqrt(1/2)*(randn(1,L)+1j*randn(1,L));
h = fftshift(h);
H = diag(h);

noise = 0.1*randn(1,J);
epsilon = 0.5;
%x = reshape(signal,num_subands,[]);%should be L long

%construct the A matrix
S = randsrc(J,L);
theta = exp(-1i*2*pi/L);
F = theta.^([0:L-1]'*[-L0:L0]);
np = 1:L0;
nn = (-L0):1:-1;
% This is for digital input only. Note that when R -> infinity,
% D then coincides with that of the paper
dn = [   (1-theta.^nn)./(1-theta.^(nn))/(L)      1/L    (1-theta.^np)./(1-theta.^(np))/(L)];
D = diag(dn);
A = S*abs(H)
A= conj(A);

y = A*s';%+ noise';

%[~,p] = chol(A*A');

x0 = A'*y;

[s history] = basis_pursuit(A,y,1.0,1.0)

figure
plot(real(s))