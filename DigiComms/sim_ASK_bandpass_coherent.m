%sim_ASK_bandpass_coherent.m 
% simulates a bandpass BASK signaling in Fig.7.2.1(coherent)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=1; M=2^b;    % # of bits per symbol and modulation order
Tb=1; Ts=b*Tb; % Bit/Symbol time 
Ns=40;         % # of sample times in Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
Es=2; A=sqrt(Es); % Energy of signal waveform
% Bandpass ASK Waveform by Eq.(7.1.1)
wc=10*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T; tt=[0:LB-1]*T;
ss=[0; 1]; su(1,:)=sqrt(2/Ts)*cos(wc*t); suT=su*T; % Eq.(7.1.2)
sw=[A*su(1,:); zeros(1,Ns)]; 
SNRdBs=[1:12]; MaxIter=10000; % Range of SNRbdB and # of iterations
for iter=1:length(SNRdBs) % Iteration for SNRdBs
   SNRbdB= SNRdBs(iter);  SNR=10^(SNRbdB/10);
   N0=2*(Es/b)/SNR;  sigma2=N0/2;  sgmsT=sqrt(sigma2/T);
   st=zeros(1,LB); y=zeros(1,LB); r=zeros(1,LB); yr=zeros(1,LB);
   nobe= 0; % Number of bit errors to be accumulated
   thrshld= sqrt(Es)/2; %+sigma)/2;
   for k=1:MaxIter
      i= ceil(rand*M); s=ss(i); % Signal index, Data bits to transmit
      for n=1:Ns  % Operation per symbol time
         st=[st(2:LB) sw(i,n)]; % Uncontaminated received signal
         yn= suT*st(LBN1:LB)'; % Correlator output with no noise
         y=[y(:,2:LB) yn]; % Correlator outputs with no noise
         wct=wcT*(n-1); 
         bp_noise= randn*cos(wct)-randn*sin(wct);
         r=[r(2:LB) sw(i,n)+sgmsT*bp_noise]; % Received signal
         yrn= suT*r(LBN1:LB)'; % Correlator output with noise
         yr=[yr(:,2:LB) yrn]; % Correlator output - DTR input
      end
      %Detector(DTR)
      if iter==10
        subplot(223), hold on
        if s==0, plot(k,yrn(1),'o');  else plot(k,yrn(1),'*'); end      
      end
      d= (yrn(1)<thrshld); % Detected data bits
      nobe = nobe+(s~=d); % Number of bit errors in a detected symbol
      if nobe>100; break; end
   end
   pobe(iter)= nobe/(k*b);
end   
SNRdBt=0:0.1:12;   SNRt=10.^(SNRdBt/10);
poe_on_theory= Q(sqrt(SNRt/4));  %Eq.(7.1.14)
subplot(224)
semilogy(SNRdBt, poe_on_theory, 'k-', SNRdBs, pobe, 'b*')
title('Error Probability for bandpass ASK Signaling')
subplot(221), plot(tt,st,'k-', tt,r,'b:')
title('Received signal r(t)')
subplot(222), plot(tt,y(1,:),'k-', tt,yr(1,:),'b:')
hold on, plot([1:4]*Tb-T,yr([1:4]*Ns),'+')
title('Output of correlator 1')