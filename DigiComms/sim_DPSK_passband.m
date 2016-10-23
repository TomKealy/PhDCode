%sim_DPSK_passband.m
% simulates a digital communication system 
%  with QDPSK signal waveforms  in Fig. 7.9
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=2; M=2^b; % # of bits per symbol and the modulation order
SNRbdBt=0:0.1:10;  SNRbt=10.^(SNRbdBt/10);
pobet= (1+(b>1))/b*Q(sqrt(b*SNRbt/2)*sin(pi/M)); %Eq.(7.4.8)
Tb=1; Ts=b*Tb; % Bit/Symbol time
Nb=16; Ns=b*Nb; % # of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size 
Es=2; sqEs=sqrt(Es); % Energy of signal waveform
% A set of b-bit symbols
ss=[0 0; 0 1; 1 1; 1 0]; phases=[0:M-1]*2*pi/M; gs='>^<vs*x+d';
wc=8*pi/Ts; t=[0:Ns-1]*T; wcT=wc*T; % Carrier frequency, etc
nd=1; % Number of delay samples
% QDPSK signal waveforms
for m=1:M, sw(m,:)= sqrt(2*Es/Ts)*cos(wc*t+phases(m));  end
su= sqrt(2/Ts)*[cos(wc*t); -sin(wc*t)]; suT= su*T;
SNRdBs=[1:10]; MaxIter=10000; % Range of SNRbdB and # of iterations
for iter=1:length(SNRdBs)
   SNRbdB= SNRdBs(iter);  SNRb=10^(SNRbdB/10);
   N0=2*(Es/b)/SNRb; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
   sws= zeros(1,LB); yr= zeros(2,LB); 
   nobe=0; % Number of bit errors to be accumulated
   is0=1; % Initial signal
   th0=1; % Initial guess of signal phase, which is possibly wrong
   for k=1:MaxIter
      i= floor(rand*M); s=ss(i+1,:); % Data bits to transmit
      is= mod(is0+i,M); is1=is+1; % Index of signal to transmit
      for n=1:Ns % Operation per symbol time
         sws=[sws(2:LB) sw(is1,n)]; 
         wct= wcT*(n-1); 
         bp_noise= randn*cos(wct)-randn*sin(wct);
         rn= sws(end-nd)+sgmsT*bp_noise; % A sample of received signal
         yr= [yr(:,2:LB) suT(:,n)*rn]; % Multiplication for correlator 
      end
      ycsk= sum(yr(:,LBN1:LB)'); % Correlator output ? DTR inut
      %Detector(DTR)
      if iter==8&k<300 
        subplot(221), hold on, plot(ycsk(1),ycsk(2),gs(i+1));
      end
      th=atan2(ycsk(2),ycsk(1));  dth=th-th0;
      if dth<-pi/M, dth=dth+2*pi; end
      [themin,mmin]=min(abs(dth-phases)); d=ss(mmin,:); % Detected bits
      nobe=nobe+sum(s~=d); if nobe>100, break; end
      is0=is; th0=th; % update the previous signal and theta
   end
   pobe(iter)= nobe/(k*b);
end   
plot(sqEs*[1 0 -1 0],sqEs*[0 1 0 -1],'ro')
subplot(222), semilogy(SNRbdBt, pobet, 'k-', SNRdBs, pobe, 'b*')
title('Probability of Bit Error for (4-ary) QDPSK Signaling')