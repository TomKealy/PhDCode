%sim_DSB_AMSC.m
% simulates DSB-AM for the message signal [4,-8,-4]
clear, clf
Ac=1; fc=50; wc=2*pi*fc;       % Amplitude/Frequency of carrier
Tb=0.1;                        % Bit interval time
T=1/fc/8; Fs=1/T;              % Sampling period/frequency
Nb=Tb/T; lt=2^(nextpow2(3*Nb)); t=[1:lt]*T; % time vector
m= ones(Nb,1)*[4 -8 -4]; m=m(:).';           % message signal m(t)
m=[m, zeros(1,lt-length(m))];
m_dsb=Ac*m.*cos(wc*t);         % AM-DSB signal (3.1.1)
y_dsb=m_dsb*2/Ac.*cos(wc*t);   % Demodulated signal (3.1.5)
N=20; Bd= fir1(N,fc*T); Ad=1;  % 20th-order FIR LPF
y_dtr=filter(Bd,Ad,y_dsb);     % Detected with LPF
plot_MOD(T,lt,m,m_dsb,y_dsb,'DSB-AMSC',y_dtr,Bd,Ad)
