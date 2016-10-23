%sim_FM.m
% simulates FM for the message signal [4,-8,-4]
clear, clf
Ac=1; fc=50; wc=2*pi*fc;            % Amplitude/Frequency of carrier
Tb=0.1;                             % Bit interval time
T=1/fc/8; Fs=1/T;                   % Sampling period/frequency
Nb=Tb/T; lt=2^(nextpow2(3*Nb)); t=[1:lt]*T; % Time vector
m= ones(Nb,1)*[4 -8 -4]; m=m(:).'; 
m=[m, zeros(1,lt-length(m))];       % Message signal m(t)
int_m(1)=0;  T2=T/2;
for i=1:lt-1
   int_m(i+1)=int_m(i)+(m(i)+m(i+1))*T2; % Integration of m(t)
end 
kF=30; % Deviation constant
m_FM=Ac*cos(wc*t+kF*int_m); % FM signal
% Demodulation/Detection
th=unwrap(angle(hilbert(m_FM)))-wc*t;  
% Eq.(1.9.8b) to get the phase of low-pass equivalent
y_FM=[0 diff(th)/T/kF];  % Demodulated/Detected signal
plot_MOD(T,lt,m,m_F,y_FM,'FM')
