%sim_FSK_bandpass_coherent.m
% simulates a digital communication system 
%  with coherent FSK signaling  in Fig.7.4.1
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=1; M=2^b; % # of bits per symbol and the modulation order
SNRbdBt=0:0.1:10;   SNRbt=10.^(SNRbdBt/10);
ss=[0; 1]; % Set of symbols
pobet= Q(sqrt(SNRbt/2)); % Eq.(7.2.9) 
Tb=1; Ts=b*Tb;           % Bit/Symbol time
Ns=40; T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
Es=2; % Energy of signal waveform
% BFSK signal waveforms
wc=8*pi; dw=2*pi; w=wc+[0:M-1]*dw; wcT=wc*T;
t=[0:Ns-1]*T; tt= [0:LB-1]*T;
nd=0; % Number of delay samples
td=0; % td=nd*T; % Delay time (time difference between XMTR and RCVR)
for m=1:M
   su(m,:)=sqrt(2/Ts)*cos(w(m)*(t-td)); 
   sw(m,:)=sqrt(2*Es/Ts)*cos(w(m)*t);
end
suT= su*T;
SNRdBs=[1:10]; MaxIter=10000; % Range of SNRbdB and # of iterations
for iter=1:length(SNRdBs) 
   SNRbdB= SNRdBs(iter);  SNR=10^(SNRbdB/10);
   N0=2*(Es/b)/SNR; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
   st=zeros(1,LB); y=zeros(M,LB); r=zeros(1,LB); yr=zeros(M,LB);
   nobe=0; % Number of bit errors to be accumulated
   for k=1:MaxIter     
      i= ceil(rand*M); s=ss(i,:); % Signal index, Data bits to transmit       
      for n=1:Ns  % Operation per symbol time  
         st=[st(2:LB) sw(i,n)]; % Uncontaminated received signal
         yn=suT*st(LBN1:LB)'; % Correlator output sample with no noise
         y=[y(:,2:LB) yn];
         wct=wcT*(n-1); 
         bp_noise= randn*cos(wct)-randn*sin(wct);
         r=[r(2:LB) sw(i,n)+sgmsT*bp_noise]; % Received signal
         yrn= suT*r(LBN1:LB)'; % Correlator output sample with noise
         yr=[yr(:,2:LB) yrn]; % noisy Correlator output buffer
      end
      % Detector(DTR)
      [yremax,lmax]= max(yrn); d=ss(lmax,:); % Detected data bits
      nobe=nobe+sum(s~=d);  if nobe>100;  break; end
   end
   pobe(iter)= nobe/(k*b);
end   
subplot(221), plot(tt,st,'k', tt,r,'b:'), title('r(t)')
subplot(222), plot(tt,y(1,:),'k', tt,yr(1,:),'b:'), hold on plot([1:4]*Tb-T,yr(1,[1:4]*Ns),'+'), title('Correlator 0 output')
subplot(224), plot(tt,y(2,:),'k', tt,yr(2,:),'b:'), hold on plot([1:4]*Tb-T,yr(2,[1:4]*Ns),'+'), title('Correlator 1 output')
subplot(223), semilogy(SNRbdBt,pobet,'k-', SNRdBs,pobe,'b*')
title('Bit Error Probability for coherent binary FSK Signaling')
