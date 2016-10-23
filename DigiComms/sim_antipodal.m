%sim_antipodal.m to simulate a binary communication system 
%  with antipodal(bipolar) signal waveforms in Fig. 5.6
clear, clf
b=1; M=2^b; % Number of bits and Number of symbols (Modulation order)
Tb=1; Ts=b*Tb;  % Bit/Symbol duration
Nb=20; Ns=b*Nb; % # of subintervals (grids) in Tb and Ts
T=Tb/Nb; % Length of subinterval 
LB=4*Ns; LBN1=LB-Ns+1; % Length of buffer in number of T
Es=2; % Energy of signal waveform
% Antipodal Waveform by Eq.(5.2.1)
su(1,:)=ones(1,Ns)/sqrt(Ts);  su(2,:)=-su(1,:); 
sw=sqrt(Es)*su; suT=su*T; ss=[0; 1];
SNRdBs= [1:10]; % Range of SNRdB
N_Iter=10000; % # of iterations for getting the probability of error
for iter=1:length(SNRdBs)
  SNRdB=SNRdBs(iter);  SNR=10^(SNRdB/10);
  sigma2=(Es/b)/SNR;  sgmsT=sqrt(sigma2/T);
  st= zeros(1,LB); y= zeros(2,LB); r= zeros(1,LB); yr= zeros(2,LB);
  nobe= 0; % Number of bit errors to be accumulated
  for k=1:N_Iter
    i=ceil(rand*M); s=ss(i); % Signal index, data bit to transmit
    for m=1:Ns
       st=[st(2:LB) sw(i,m)]; % Uncontaminated received signal
       y=[y(:,2:LB) suT*st(LBN1:LB)']; % No noise
       r=[r(2:LB) sw(i,m)+sgmsT*randn]; % Received signal
       yr=[yr(:,2:LB) suT*r(LBN1:LB)']; % Correlator output-DTR input
    end
    % Detector(DTR)
    d=(yr(1,end)<0); % Detected data bit
    % if yr(1,end)>0, d=0; else d=1; end
    nobe = nobe + sum(d~=s); if nobe>100, break; end
  end
  pobe(iter)=nobe/(k*b); % probability of bit detection errors
end   
SNRdBt=0:0.1:10;  SNRt=10.^(SNRdBt/10);
pobe_on_theory= Q(sqrt(SNRt));  % Eq.(5.2.34a)
subplot(223)
semilogy(SNRdBt, pobe_on_theory, 'k-', SNRdBs, pobe, 'b*')
set(gca,'fontsize',9), title('Probability of Error for Antipodal Signaling')
t= [0:LB-1]*T; % time duration for buffering
subplot(221), plot(t,st,'k:', t,r,'b-')
set(gca,'fontsize',9), title('Received signal r(t)')
subplot(222), plot(t,y(1,:),'k:', t,yr(1,:),'b-')
set(gca,'fontsize',9), title('Output of correlator 1')
subplot(224), plot(t,y(2,:),'k:', t,yr(2,:),'b-')
set(gca,'fontsize',9), title('Output of correlator 2')
