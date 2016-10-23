%sim_FSK_bandpass_noncoherent.m 
% simulates a digital communication system 
%  with non-coherent FSK signaling  in Fig.7.4.2
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=1; M=2^b; % # of bits per symbol and the modulation order
SNRbdBt=0:0.1:10; SNRbt=10.^(SNRbdBt/10); 
pobet=exp(-SNRbt/4)/2;  % Eq.(7.2.20) 
Tb=1; Ts=b*Tb;          % Bit/Symbol time
Ns=32; T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
Es=2; % Energy of signal waveform
wc=8*pi/Ts; dw=2*pi/Ts; w=wc+[0:M-1]*dw; wcT=wc*T;
t=[0:Ns-1]*T;  tt=[0:LB-1]*T;
nd=0; % Delay sample time (time offset between XMTR and RCVR clocks)
for m=1:M  % BFSK signal waveforms
   sw(m,:)=sqrt(2*Es/Ts)*cos(w(m)*t);
   su(2*m-1:2*m,:)=sqrt(2/Ts)*[cos(w(m)*t); sin(w(m)*t)];
end
suT=su*T; ss=[0; 1];
SNRdBs=[1:10]; MaxIter=10000; % Range of SNRbdB and # of iterations
for iter=1:length(SNRdBs)
SNRbdB= SNRdBs(iter);  SNR=10^(SNRbdB/10);
N0=2*(Es/b)/SNR; sigma2=N0/2; sgmsT=sqrt(sigma2/T);   
st=zeros(1,LB); sws=zeros(1,LB);  y=zeros(M,LB); yn=zeros(M,1);
r=zeros(1,LB);  yr=zeros(M,LB);   yrn=zeros(M,1);
nobe= 0; % Number of bit errors to be accumulated
for k=1:MaxIter
i=ceil(rand*M); s=ss(i,:); % Signal index, Data bits to transmit   for n=1:Ns  % Operation per symbol time
sws=[sws(2:LB) sw(i,n)]; % updates the signal waveform buffer
st=[st(2:LB) sws(end-nd)]; % Uncontaminated received signal
ycss=suT*st(LBN1:LB)'; % Correlator output with no noise
for m=1:M,  yn(m)=ycss(2*m-1)^2+ycss(2*m)^2;  end
y=[y(:,2:LB) yn]; % Envelope detector output
wct=wcT*(n-1); bp_noise= randn*cos(wct)-randn*sin(wct);
r=[r(2:LB) sws(end-nd)+sgmsT*bp_noise]; % Received signal
ycss=suT*r(LBN1:LB)'; % Correlator output with noise
for m=1:M,  yrn(m)=ycss(2*m-1)^2+ycss(2*m)^2;  end
yr=[yr(:,2:LB) yrn]; % Envelope detector output buffer
end
     %Detector(DTR)
     [yremax,mmax]=max(yrn);  d=ss(mmax,:); % Detected data bits
   nobe=nobe+sum(s~=d);  if nobe>100;  break; end
end
pobe(iter)= nobe/(k*b);
end   
subplot(221), plot(tt,st,'k', tt,r,'b:'), title(' r(t)')
subplot(222), plot(tt,y(1,:),'k', tt,yr(1,:),'b:')
hold on, plot([1:4]*Tb-T,yr(1,[1:4]*Ns),'+'), title('y0(t)')
subplot(224), plot(tt,y(2,:),'k', tt,yr(2,:),'b:')
hold on, plot([1:4]*Tb-T,yr(2,[1:4]*Ns),'+'), title('y1(t)')
subplot(223), semilogy(SNRbdBt,pobet,'k-', SNRdBs,pobe,'b*')
title('Bit Error Probability for non-coherent BFSK Signaling')
