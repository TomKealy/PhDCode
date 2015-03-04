clear all
close all

% W=6000 bin singal divided into L=500 subbands of bandwidth B=12 bins
% (n.b. all in Fourier domain).


%number of samples per period
s=4;

%RF frequency
f=4e9;

%pulse repetition frequency
prf=1/30e-9;

%sampling frequency
fs=s*f;

%Total Simulation time
T=30e-9;

t=0:1/fs:T;

%generating pulse train
x=pulstran(t,15e-9,'gauspuls',f,0.5);

%length of the signal
W=length(x);

L=239; %number of sub-bands
L0 = (L-1)/2; %re-labelling via shift
J = 50; %number of nodes
B=W/L; %bandwidth per subband
K = 0.02*L; % sparsity

figure;
subplot(2,1,1);
plot(t,x)
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Original Signal, UWB Pulse RF freq=%g GHz',f/1e9));

%taking Discrete time Fourier Transform of the signal
xf=fft(x);

xfmag=10*log10(abs(xf));

subplot(2,1,2);
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Discrete Fourier Transform of UWB pulse');

noise = 0.1*randn(1,J);
epsilon = 0.5;

%construct the A matrix

h = sqrt(1/2)*(randn(1,L)+1j*randn(1,L));
h = fftshift(h);
H = diag(h);
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
A= conj(A);

%taking random frequency measurements
y=(A*x');

% Calculating Initial guess
x0=A'*y;

%Running the recovery Algorithm
cvx_begin
    variable x(N) complex;
    minimise (norm(xf,1))
    subject to
        A*xf' == y
cvx_end

s = ifft(s)
figure;
subplot(2,1,1)
plot(t,x)
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Original Signal, UWB Pulse RF freq=%g GHz',f/1e9));

subplot(2,1,2)
plot(t,real(xf),'r')
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Recovered UWB Pulse Signal with %d random samples',K));