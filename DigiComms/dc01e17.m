%dc01e17.m
% Sinusoidal-modulated FM signal
clear, clf
beta=0:.05:15;
[J15,JJ]=Jkb(15,beta);
for k=0:15
   plot(beta,JJ(k+1,:)), hold on
end 
title('Bessel ftn of the 1st kind & kth order')
A=1; beta=5; 
wc=2000*pi; wm=40*pi; % center frequency and swaying frequency 
T=0.0002; N=500; t=[0:N-1]*T;  
st=A*cos(wc*t+beta*sin(wm*t));
tt=[0:400]*0.00005; stt=A*cos(wc*tt+beta*sin(wm*tt));
pause, clf, subplot(311), plot(tt,stt)
Sk=fftshift(fft(st)); Sk=[Sk Sk(1)]*T; %DFT spectrum made symmetric
subplot(312), plot([-N/2: N/2], abs(Sk))
P=N*T; w0=2*pi/P;  Kc=wc/w0; M=wm/w0;
kk=[Kc-20:Kc+20]; %the band around the center frequency
subplot(313), stem(kk, abs(Sk(kk+N/2+1))) %the LHS of Eq.(E1.17.6)
hold on, pause
Jk= zeros(size(Sk));
for k=0:10
   Jk(N/2+k*M+Kc+1)=Jkb(k,beta); Jk(N/2+k*M-Kc)=Jkb(k,beta);
   if k>0
     Jk(N/2-k*M+Kc+1)=(-1)^mod(k,2)*Jk(N/2+k*M+Kc+1);
   Jk(N/2-k*M-Kc)=(-1)^mod(k,2)*Jk(N/2+k*M-Kc);
 end
end
stem(kk, P/2*A*abs(Jk(kk+N/2+1)),'r') %the RHS of Eq.(E1.17.6)
