%sim_PSK_passband.m
% simulates a digital communication system 
%  with QPSK signal waveforms in Fig.7.7
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=2; M=2^b; % # of bits per symbol and the modulation order
SNRbdBt=0:0.1:10;   SNRbt=10.^(SNRbdBt/10);
%pobet=Q(sqrt(SNRbt)).*(1-0.5*Q(sqrt(SNRbt))); %Eq.(5.2.59) in [P-4]
pobet= (1+(b>1))/b*Q(sqrt(b*SNRbt)*sin(pi/M)); %Eq.(7.3.7)
Tb=1; Ts=b*Tb; % Bit/Symbol times 
Nb=32; Ns=b*Nb; % # of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size 
Es=2; sqEs=sqrt(Es); % Energy of signal waveform
% A set of b-bit symbols
ss=[0 0; 0 1; 1 1; 1 0]; 
%g_code=gray_code(b); for m=1:M, ss(m,:)=deci2bin(g_code(m),b); end
phases=2*pi/M*[0:M-1];  gs='>^<vs*x+d';
% Carrier frequency
wc=8*pi/Ts; t=[0:Ns-1]*T; wcT=wc*T; 
nd=0; % Number of delay samples 
% QPSK signal waveforms
for m=1:M, sw(m,:)= sqrt(2*Es/Ts)*cos(wc*t+phases(m));  end
su= sqrt(2/Ts)*[cos(wc*t); -sin(wc*t)]; suT= su*T;
SNRdBs=[1:10]; MaxIter=10000; % Range of SNRbdB and # of iterations 
for iter=1:length(SNRdBs)
   SNRbdB= SNRdBs(iter);  SNR=10^(SNRbdB/10);
   N0=2*(Es/b)/SNR;  sigma2=N0/2;  sgmsT=sqrt(sigma2/T);
   sws= zeros(1,LB);  yr= zeros(2,LB);
   nobe= 0; % Number of bit errors to be accumulated
   for k=1:MaxIter
      i=ceil(rand*M); s=ss(i,:); % Signal index, Data bits to transmit       
      for n=1:Ns % Operation per symbol time
         sws=[sws(2:LB) sw(i,n)]; % Buffering the signal waveform
         wct= wcT*(n-1); bp_noise= randn*cos(wct)-randn*sin(wct);
         rn= sws(end-nd)+sgmsT*bp_noise; % A sample of received signal
         yr= [yr(:,2:LB) suT(:,n)*rn]; % Multiplication for correlator
      end
      ycsk= sum(yr(:,LBN1:LB)')'; % Correlator output - DTR input
      %Detector(DTR)
      if iter==9&k<300  % Signal constellation diagram
        subplot(221), hold on, plot(ycsk(1),ycsk(2),gs(i))
      end
      th= atan2(ycsk(2),ycsk(1));
      if th<-pi/M, th=th+2*pi; end
      [themin,lmin]= min(abs(th-phases)); % Eq.(7.3.9)
      d= ss(lmin,:); % Detected data bits
      nobe=nobe+sum(s~=d); % Number of bit errors in a detected symbol
      if nobe>100, break; end
   end
   pobe(iter)= nobe/(k*b);
end   
plot(sqEs*[1 0 -1 0],sqEs*[0 1 0 -1],'ro')
subplot(222), semilogy(SNRbdBt,pobet,'k-', SNRdBs,pobe,'b*')
title('Probability of Bit Error for (4-ary) QPSK Signaling')