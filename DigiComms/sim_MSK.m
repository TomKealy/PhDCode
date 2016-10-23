%do_MSK.m 
% simulates a digital communication system with MSK signaling
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
b=1; M=2^b; % # of bits per symbol and modulation order
SNRbdBt=0:0.1:12;   SNRbt=10.^(SNRbdBt/10);
pbe_PSK= prob_error(SNRbdBt,'PSK',b,'bit');
pbe_FSK= prob_error(SNRbdBt,'FSK',b,'bit');
Tb=1; Ts=b*Tb; % Bit/Symbol time 
Nb=16; Ns=Nb*b; % # of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN2=LB-2*Ns+1; % Sample time and Buffer size
Es=2;  % Energy of signal waveform
% QPSK signal waveforms
ds=[0 1];
wc=2*pi; wcT=wc*T; 
t=[0:2*Nb-1]*T; pihTbt=pi/2/Tb*t;
su=sqrt(2/Tb)*[cos(pihTbt-pi/2).*cos(wc*(t-Tb)); ??????????????????????]; % Eq.(P7.6.8)
su=[fftshift(su(1,:)); su(2,:)];  suT=su*T;
ik=[1 2 1 2 1 1 2 2 2 1];
sq2EbTb=sqrt(2*Es/b/Tb); pihTbT=pi/2/Tb*T;
SNRdBs=[1 3 10]; % Range of SNRbdB
MaxIter=10000; % Number of iterations
for iter=1:length(SNRdBs) 
   SNRbdB=SNRdBs(iter);  SNR=10^(SNRbdB/10);
   N0=2*(Es/b)/SNR; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
   yr=zeros(2,LB); yc=0; ys=0; th0=0; t=0; wct=0; tht=th0; mn=1;
   nobe= 0; % Number of bit errors to be accumulated
   for k=1:MaxIter
      i=ceil(rand*M); s(k)=ds(i); sgn=s(k)*2-1; 
      for m=1:Nb % Operation per bit time
         bp_noise= randn*cos(wct)-randn*sin(wct);
         rn= sq2EbTb*cos(wct+tht) + sgmsT*bp_noise;
         yr=[yr(:,2:LB) suT(:,mn)*rn]; % Correlator output - DTR input
         t=t+T; wct=wct+wcT; tht=tht+sgn*pihTbT; mn=mn+1; 
      end
      if tht<-pi2, tht=tht+pi2; elseif tht>pi2, tht=tht-pi2; end
      thtk(k+1)=tht; 
      if mn>size(su,2), mn=1; end
      if mod(k,2)==1
        yc=sum(yr(?,LBN2:LB));  th_hat(k)=??*(yc<0);
        % if k>1&th_hat(k-1)<0, th_hat(k)=-th_hat(k); end
       else 
        ys=sum(yr(?,LBN2:LB));  th_hat(k)=????*(2*(ys<0)-1); 
      end
      if k>1 % Detector(DTR)
        d=(dth_hat>0);  if abs(dth_hat)>=pi, d=~d; end
        nobe = nobe+(d~=s(???)); if nobe>100; break; end
      end         
   end
   pobe(iter)= nobe/((k-1)*b);
end   
semilogy(SNRbdBt,pbe_PSK,'k', SNRbdBt,pbe_FSK,'k:', SNRdBs,pobe,'b*')