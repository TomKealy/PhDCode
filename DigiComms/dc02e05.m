%dc02e05.m: Statistical Property of White Noises
% the autocorrelation and distribution (histogram) of a white noise
clear, clf
Nb=8; %number of samples per bit time interval
M=5000; N=M*Nb; % Sample size
Tb=1; dt=Tb/Nb; % Bit time interval and sampling interval
xlim=10; xx=-xlim:0.02:xlim; 
NB=40; % Number of bins
randn('state',0); %return randn() to its default initial state 
w=randn(1,N)/sqrt(dt); % divide by sqrt(dt) so that integral(w^2(t))dt=1
for i=1:2
   if i==1
     sigma2= w*w'/N % phiw(0)=integral(w^2(t))dt/(N*dt)
     sgm2=1/dt; f=exp(-xx.^2/2/sgm2)/sqrt(2*pi*sgm2);
     nn=[-N/2:N/2]; tt=nn*dt; nnt=nn+N+1;
    else
     w=cos(2*pi*[0:Nb-1]/Nb)*reshape(w,Nb,M)*sqrt(2/Tb)*dt; %Eq.(E2.5.1)
     sigma2= w*w'/M % phiw[0]=sum(w^2(n))*Tb/(M*Tb)    
     f=exp(-xx.^2/2)/sqrt(2*pi);
     nn=[-M/2:M/2]; tt=nn*Tb; nnt=nn+M+1;
   end
   phiw=xcorr(w,'u');
   subplot(219+2*i), plot(tt,phiw(nnt))
   idx=find(abs(w)<=xlim); 
   w1=w(idx); N1=length(w1);
   subplot(220+2*i)
   [ns,cs]=hist(w1,NB); 
   dx=cs(2)-cs(1);
   hist(w1,NB), hold on, plot(xx,N1*f*dx,'r')
   % plot(xx,f*dx,'r',cs,ns/N1,':')
end
