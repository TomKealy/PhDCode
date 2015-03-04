%Priamry User Signal
path(path,'./Optimization')

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

s = zeros(L,1); %psd
s(positions) = 1;
sprime = s;
figure
plot(s)

%Secondary User system

h = sqrt(1/2)*(randn(1,L)+1j*randn(1,L));
h = fftshift(h);
H = diag(h);

noise = 0.001*randn(1,J);

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
A = S*F*D*H;
A = conj(A);

y = A*s+ noise';

cvx_begin
    variable s(L) complex;
    minimise (norm(s,1))
    subject to
        A*s == y
cvx_end


[solution] = GPSR_BB(y, A, 0.3, 'Verbose', 0, 'StopCriterion', 3, 'ToleranceA', 1e-6);

opts = spgSetParms('verbosity',0);
d = spg_bp(A, y, opts);

figure
plot(20*log10(abs((d))))

figure
plot(20*log10(abs((solution))))

figure
plot(20*log10(abs(s)))

mse_spg = norm(d - sprime)/norm(s)
mse_gpsr = norm(solution - sprime)/norm(s)
mse_cvx = norm(s-sprime)/norm(s)