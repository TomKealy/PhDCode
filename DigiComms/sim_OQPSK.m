%do_OQPSK.m 
% simulates a digital communication system 
%  with O(ffset)QPSK signal waveforms in Fig.P7.4
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=2; M=2^b;
SNRbdBt=0:0.1:10;  SNRbt=10.^(SNRbdBt/10);
pobet=prob_error(SNRbdBt,'PSK',b,'bit');
Tb=1; Ts=b*Tb; % Bit/symbol time
Nb=16; Ns=Nb*b; % # of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
Es=2; sqEb=sqrt(Es/b); % Energy of signal waveform
% QPSK signal waveforms
ss=[0 0; 0 1; 1 1; 1 0];
wc=2*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T; 
su=sqrt(2/Ts)*???????????????????????;  suT=su*T;
sw=sqEb*su;
SNRdBs=[1:3:10]; MaxIter=10000; % Range of SNRdB, # of iterations
for iter=1:length(SNRdBs)
   SNRbdB=SNRdBs(iter);  SNRb=10^(SNRbdB/10);
   N0=2*(Es/b)/SNRb; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
   yr= zeros(2,LB); yc=zeros(1,2); ys=zeros(1,2);
   iq=0; % Initialize the quadrature bit arbitrarily
   nobe=0; % Number of bit errors to be accumulated
   for k=1:MaxIter
      i=ceil(rand*M); s=ss(i,:);  wct=-wcT; 
      for n=1:b % Operation per symbol time
         if n==1, ii=2*ss(i,1)-1; % In-phase bit
          else iq=2*ss(i,2)-1; % Quadrature bit
         end
         mn=0;
         for m=1:Nb % Operation per bit time
            wct= wct+wcT; mn=mn+1;
            bp_noise= randn*cos(wct)-randn*sin(wct);
            rn=ii*sw(1,mn)+iq*sw(2,mn)+sgmsT*bp_noise;
            yr=[yr(:,2:LB) suT(:,mn)*rn]; % Multiplier
         end    
         if n==2  % sampled at t=2k*Tb
           yc=[yc(2) sum(yr(?,LBN1:LB))]; % Correlator output
          else  % sampled at t=(2k-1)*Tb
           ys=[ys(2) sum(yr(?,LBN1:LB))]; % Correlator output
         end
      end
      d=([yc(?) ys(?)]>0); % Detector(DTR)
      if k>1, nobe=nobe+sum(s0~=d); end
      if nobe>100, break; end
      s0= s;
   end
   pobe(iter)= nobe/(k*b);
end   
pobe
semilogy(SNRbdBt,pobet,'k-', SNRdBs,pobe,'b*')
title('Probability of Bit Error for (4-ary) QPSK Signaling')
