%sim_PM.m
% simulates PM for the message signal [4,-8,-4]
clear, clf
Ac=1; fc=50;  wc=2*pi*fc; % Amplitude/Frequency of carrier
Tb=0.1; ts=1/fc/8; fs=1/ts; %bit interval, sampling period/freq
Nb=Tb/ts; lt=2^(nextpow2(3*Tb/ts)); t=[1:lt]*ts; % time vector
m=4*[ones(1,Nb), -2*ones(1,Nb), -ones(1,Nb)]; % Message signal
m=[m, zeros(1,lt-length(m))];  
kP=0.01;%deviation constant 30*ts
kP=30;
m_PM=Ac*cos(wc*t+kP*m); %PM signal
%demodulation
th=unwrap(angle(hilbert(m_PM)))-wc*t; %phase of analytic signal
y_PM=th/kP; % demodulated signal
plot_MOD(ts,lt,m,m_PM,y_PM,'PM')

t=[1:lt]*ts;
m_PM2=modulate(m,fc,fs,'pm',kP);
subplot(423), hold on, plot(t,m_PM2,'r')
y_PM2=demod(m_PM2,fc,fs,'pm',kP);
subplot(425), hold on, plot(t,y_PM2,'r')