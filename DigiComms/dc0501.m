%dc0501.m
% simulates binary communication systems in Figs.5.2/5.3(5.4)
%  with antipodal (bipolar) signal waveforms in Fig.5.6
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
Tb=1; sqTb=sqrt(Tb); % Tb=Bit time and its square root
Nb=20; % Number of samp[le times/grids in Tb for integration
T=Tb/Nb; % Sample time
L=6; LB=L*Nb; LBN1=LB-Nb+1; % Length of buffer in number of T
Es=2; sqEs=sqrt(Es); % Energy of signal waveform
su=[1;-1]/sqTb*ones(1,Nb); % Unit signal waveform Eq.(5.1.5)
t= [0:LB-1]*T; % Time duration for buffering
RC=Tb/1.2564; g_RC=exp(-t/RC)/RC; % Impulse response of RC filter
gr_RC=fliplr(g_RC); % Reversed/delayed impulse response of RC filter
sw=sqEs*su;   % Antipodal waveform by Eq.(5.2.1)
SNRdB=6;  SNR=10^(SNRdB/10);  N0=2*Es/SNR;  % Noise power Eq.(5.1.25)
sigma2=N0/2; sigma_sT=sqrt(sigma2/T); % Noise power (variance)
st= zeros(1,LB); y= zeros(2,LB);  y_RC= zeros(1,LB);
r= zeros(1,LB); yr= zeros(2,LB); yr_RC= zeros(1,LB);
MaxIter=20; % # of iterations 
for n=1:MaxIter
   s(n)=floor(2*rand); i=s(n)+1;
   for m=1:Nb % Bit processing
      st=[st(2:LB) sw(i,m)]; % Uncontaminated received signal
      y_RC=[y_RC(:,2:LB) gr_RC*st'*T]; % RC filter output with no noise
      y=[y(:,2:LB) su*st(LBN1:LB)'*T];% Correlator output with no noise
      r=[r(2:LB) sw(i,m)+sigma_sT*randn]; % Received signal with noise
      yr_RC=[yr_RC(:,2:LB) gr_RC*r'*T];  % RC filter output with noise
      yr=[yr(:,2:LB) su*r(LBN1:LB)'*T]; % Correlator output with noise
   end
   %Detector(DTR)
   D_RC(n)= (yr_RC(1,end)<0);
   D_c(n)= (yr(1,end)<0);
end
poe_RC=sum(s~=D_RC)/MaxIter, poe_c=sum(s~=D_c)/MaxIter
subplot(221), plot(t,st,'k:', t,r,'b-'), title('Received signal r(t)')
subplot(222), plot(t,y(1,:),'k:', t,yr(1,:),'b-')
title('Output of correlator 1')
subplot(223), plot(t,y_RC,'k:', t,yr_RC,'b-')
title('Output of RC filter')
subplot(224), plot(t,y(2,:),'k:', t,yr(2,:),'b-')
title('Output of correlator 2')
