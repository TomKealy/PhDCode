%sim_Delta_Sigma.m
clear, clf
T=0.01; nmax=200; t=[0:nmax-1]*T;
Delta2=1.5; % Half of the 1-bit quantization step size
% Design an elliptic filter
fp=0.1; fs=0.15; Rp=3; As=40; ; 
[N,fc]=ellipord(fp,fs,Rp,As); [B,A]=ellip(N,Rp,As,fc);
wi=zeros(max(length(A),length(B))-1,1);
din=0; dqn=0; yqn=0; % dn=0;
for n=1:nmax
   x(n)= sin(6*(n-1)*T); % +0.2*(rand-0.5);
   d(n)= x(n) - dqn; % Comparator input
   di(n)= din + d(n)*?; % Integrator output: di(n)=di(n-1)+d(n)*T
    % Transfer function of this integrator: Di[z]/D[z]=T/(z-1)
   dq(n)= sign(di(n))*??????;
   din=di(n);  dqn=dq(n); % Store/Delay
   [yq(n),wi]=filter(B,A,dq(n),wi);
end
subplot(611), plot(t,x,'k-'), hold on, stairs(t,di,'b:'), stairs(t,d,'b-'), stairs(t,dq,'r-'), stairs(t,yq,'m'), 
subplot(612), plot(t,x-yq,'r-')
