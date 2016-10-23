%dc01p02.m 
% Envelope and phase for lowpass equivalent of a bandpass signal
clear, clf
ts=0.0001; fs=1/ts; % Sampling Period/Frequency
N=1024; t=[-(N-1):N]*ts; % duration of signal
n1=[N-400:N+400]; t1=t(n1); % time vector for display duration
f1=200; w1=2*pi*f1; B=200; % frequency/bandwidth of signal x(t)
x= sinc(B/pi*t).*cos(w1*t+pi/6);  x1=x(n1);  
xh=imag(hilbert(x));  % imag(x(t)+j x^(t))=x^(t)
xc= x.*cos(w1*t)+xh.*sin(w1*t); 
xs=-x.*sin(w1*t)+xh.*cos(w1*t); 
env1=????????????????; % Eq.(P1.1.1a1)
env2=?????????????????; % Eq.(P1.1.1a2)
ph1=????????????; % Eq.(P1.1.1b1)
ph2=????????????????; % Eq.(P1.1.1b2)
%ph2=mod(ph2,2*pi);
%idx=find(ph2>pi);
%ph2(idx)=ph2(idx)-2*pi;
subplot(221)
plot(t1,x(n1),'k:', t1,env1(n1),'b', t1,env2(n1),'r') 
subplot(222)
plot(t1,ph1(n1),'b', t1,ph2(n1),'r')
