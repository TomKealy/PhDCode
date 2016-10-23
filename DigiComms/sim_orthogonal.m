%sim_orthogonal.m
% simulates a binary communication system 
%  with orthogonal signal waveforms in Fig.5.9
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=1; M=2^b; % # of bits/symbols
Tb=1; Ts=b*Tb;  % Bit/Symbol duration
Nb=20; Ns=b*Nb; % # of subintervals (grids) in Tb and Ts (for integration)
T=Tb/Nb; % Ssubinterval length
LB=4*Ns; LBN1=LB-Ns+1; % Length of buffer in number of T
Es=2; % Energy of signal waveform
% Orthogonal Waveform su0&s0, su1&s1 satisfying Eq.(5.2.21)
su=[1 1; 1 -1]; su=ones(Ns/M,1)*reshape(su',1,M*M);
su=reshape(su(:),Ns,M)'/sqrt(Ts);
suT=su*T;  sw=sqrt(Es)*su; 
for m=1:M, ss(m,:)=deci2bin(m-1,b); end  % ss=[0; 1];
SNRdBs= [1:10]; % Range of SNRdB
N_Iter=6000; % # of iterations for getting the error probability
for iter=1:length(SNRdBs) 
  SNRdB= SNRdBs(iter);  SNR= 10^(SNRdB/10);
  sigma2=(Es/b)/SNR;  sgmsT=sqrt(sigma2/T);
  st= zeros(1,LB);  y= zeros(M,LB);
  r= zeros(1,LB);   yr= zeros(M,LB);
  nobe= 0; % Number of errors to be accumulated
  for k=1:N_Iter
    i= ceil(rand*M); s=ss(i,:); % Signal index, data bit to transmit
    for m=1:Ns
       st=[st(2:LB) sw(i,m)]; % Uncontaminated received signal
       y=[y(:,2:LB) suT*st(LBN1:LB)'];
       r=[r(2:LB) sw(i,m)+sgmsT*randn]; % Received signal
       yre= suT*r(LBN1:LB).'; 
       yr=[yr(:,2:LB) yre]; % Correlator output - DTR input
    end
    % Detector(DTR)
    [ymax,imax]=max(yre); d = ss(imax,:); % Detected data bits
    nobe = nobe+sum(d~=s);  
    if nobe>100, break; end
  end
  pobe(iter)= nobe/(k*b);
end   
SNRdBt=0:.1:10;   SNRt=10.^(SNRdBt/10);
pobe_on_theory= Q(sqrt(SNRt/2));  % Eq.(5.2.34c)
subplot(223)
semilogy(SNRdBt,pobe_on_theory,'k-', SNRdBs,pobe,'b*')
set(gca,'fontsize',9), title('Probability of Error for Orthogonal Signaling')
t= [0:LB-1]*T; % Time duration for buffering
subplot(221), plot(t,st,'k:', t,r,'b-')
set(gca,'fontsize',9), title('Received signal r(t)')
subplot(222), plot(t,y(1,:),'k:', t,yr(1,:),'b-')
set(gca,'fontsize',9), title('Output of correlator 1')
subplot(224), plot(t,y(2,:),'k:', t,yr(2,:),'b-')
set(gca,'fontsize',9), title('Output of correlator 2')
