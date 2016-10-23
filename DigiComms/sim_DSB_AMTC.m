%sim_DSB_AMTC.m
% simulates conventional AM(DSB-AMTC) for the message signal [4,-8,-4]
clear, clf
Ac=1; fc=50; wc=2*pi*fc;      % Amplitude/Frequency of carrier
Tb=0.1;                       % Bit interval time
T=1/fc/8; Fs=1/T;             % Sampling period/frequency
Nb=Tb/T; lt=2^(nextpow2(3*Nb)); t=[1:lt]*T; % Time vector
m= ones(Nb,1)*[4 -8 -4]; m=m(:).'; % Message signal m(t)
m=[m, zeros(1,lt-length(m))];  m_max=max(abs(m));
a=0.9;                             % Demodulation parameter
m_cv=Ac*(1+a*m/m_max).*cos(wc*t);  % CV-AM signal
y_dtr=abs(hilbert(m_cv));          % Envelope detector output
y_cv=(y_dtr/Ac-1)*m_max/a;         % Demodulated 
plot_MOD(T,lt,m,m_cv,y_dtr,'DSB-AMTC',y_cv)
