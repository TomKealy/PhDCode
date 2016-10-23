%sim_DPCM.m
% adaptive predictor in XMTR and the predictor 
%  with the same (transmitted) coefficients in RCVR
clear, clf
T=0.01; nmax=200; t=[0:nmax-1]*T;  
%b0=min(x); bN=max(x); xmax=max(abs([b0 bN]))
b0=-0.4; bN=0.4; N=8; %left/rightmost boundary  & # of bits/sample
delta=(bN-b0)/N;  %the width of the quantization interval
b=b0+[0:N]*delta; %the boundary vector
c=b(1:N)+delta/2; %quantization level vector
M=2;  M1=M+1; % dimension of predictor(RLSE)
% initialize the parameters of the predictor
P=eye(M1,M1); a=zeros(M1,1); alpha=1;  
%initialize the (memory) buffers of XMTR/predictor
xqs=[1; zeros(M,1)];  xqhn=0; dn=0; dqn=0;
%initialize the (memory) buffers of RCVR/predictor
yqs=[1; zeros(M,1)]; 
for n=1:nmax
   x(n)= sin(6*(n-1)*T);
   xqs= [1; xqs(3:M1); xqhn+dqn]; %the predictor output history
   if n<40 %training period for adaptive predictor
     [xqh(n),P,a]=prdctr(dn,xqs,M,P,a,alpha); %adaptive predictor
    else    
      xqh(n)= a'*xqs;  %non-adaptive predictor
   end
   d(n)=x(n)-xqh(n); %prediction error
   dq(n)=adc(d(n),b,c); %quantizer at XMTR(encoder)
   dn=d(n); dqn=dq(n); xqhn=xqh(n); %store/delay
   yq(n)=a'*yqs +dqn; %the output of predictor at RCVR(decoder)
   yqs=[1; yqs(3:M1); yq(n)];  %store/delay
end
subplot(211), plot(t,x,'k-'), hold on
stairs(t,d,'b-'),  stairs(t,dq,'r-'), stairs(t,yq,'m-')
stairs(t,d,'b-'),  stairs(t,dq,'r-'), stairs(t,yq,'k-')
axis([0 2 -1.15 1.15]), set(gca,'fontsize',9)
subplot(212), plot(t,x-yq,'r-')
axis([0 2 -0.13 0.13]), set(gca,'fontsize',9)
