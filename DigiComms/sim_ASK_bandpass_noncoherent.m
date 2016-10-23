%sim_ASK_bandpass_noncoherent.m 
% simulates a binary communication system 
%  with bandpass ASK signal waveforms in Fig.7.2.2,
%  where the envelope detector is used for non-coherent detection.
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=1; M=2^b;    % # of bits per symbol and modulation order
Tb=1; Ts=b*Tb; % Bit/Symbol time
Ns=40;         % # of sample times in Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
Es=2; A=sqrt(Es); % Energy of signal waveform
% Bandpass ASK Waveform by Eq.(7.1.1)
wc=10*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T;  wct=wc*t;  tt=[0:LB-1]*T;
ss=[0;1]; su=sqrt(2/Ts)*[cos(wct); -sin(wct)]; suT=su*T; % Eq.(7.1.2)
sw(1,:)= A*su(1,:); sw(2,:)= zeros(1,Ns);
SNRdBs=[1:12]; MaxIter=10000; % Range of SNRbdB and # of iterations  
for iter=1:length(SNRdBs)
   SNRbdB= SNRdBs(iter);  SNR=10^(SNRbdB/10);
   N0=2*(Es/b)/SNR;  sigma2=N0/2;  sgmsT=sqrt(sigma2/T);
   st=zeros(1,LB);  y=zeros(1,LB);  r=zeros(1,LB);  yr=zeros(1,LB);
   nobe= 0; % Number of errors
   thrshld= Es/2; %+sigma)/2;
   for k=1:MaxIter
      i=ceil(rand*M); s=ss(i); % Signal index, Data bits to transmit
      for n=1:Ns  % Operation per symbol time
         st=[st(2:LB) sw(i,n)]; % Uncontaminated received signal
         yn=su*st(LBN1:LB)'*T; % Quadrature correlator outputs
         y=[y(:,2:LB) yn(1)^2+yn(2)^2]; % Square device output
         wct= wcT*(n-1);
         bp_noise= randn*cos(wct)-randn*sin(wct); % Bandpass noise
         r=[r(2:LB) sw(i,n)+sgmsT*bp_noise]; % Received signal
         yrn=suT*r(LBN1:LB)'; % noisy Quadrature correlator outputs yr=[yr(:,2:LB) yrn(1)^2+yrn(2)^2]; % Envelope Detector output
      end
      % Detector(DTR)
      if iter==10
        subplot(223)
        if s==0, plot(k,yr(LB),'o');  else plot(k,yr(LB),'*');  end
        hold on
      end
      d= (yr(LB)<thrshld); % Detected data bit
      nobe = nobe+(s~=d);  if nobe>100;  break; end
   end
   pobe(iter)= nobe/(k*b);
end
SNRdBt=0:0.1:12;   SNRt=10.^(SNRdBt/10);
poe_on_theory= exp(-SNRt/8)/2; %+Q(A/2)/2;  %Eq.(7.1.21) 
subplot(224), semilogy(SNRdBt,poe_on_theory,'k-', SNRdBs,pobe,'b*')
title('Error probability for bandpass ASK Signaling')
subplot(221), plot(tt,st,'k-', tt,r,'b:'), title('r(t)')
subplot(222), plot(tt,y(1,:),'r-', tt,yr(1,:),'b:')
hold on, plot([1:4]*Tb-T,yr([1:4]*Ns),'+')
title('Envelope detector output')