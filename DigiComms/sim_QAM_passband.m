%sim_QAM_passband.m
% simulates a digital communication system in Fig.7.13
%  with QAM signal waveforms in Fig.7.11
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
b=4; M=2^b; L=2^(b/2); % # of bits per symbol and the modulation order
SNRbdBt=0:0.1:15;   SNRbt=10.^(SNRbdBt/10);
Pm=2*(1-1/L)*Q(sqrt(3/2*b*SNRbt/(M-1))); % Eq.(7.1.5)
pobet= (1-(1-Pm).^2)/b; % Eq.(7.5.8) with (7.5.6)
Tb=1; Ts=b*Tb;  % Bit/Symbol time 
Nb=16; Ns=b*Nb; % # of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
ssc=[0 0; 0 1; 1 1; 1 0]; sss=ssc; 
wc=8*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T; 
su=sqrt(2/Ts)*[cos(wc*t); -sin(wc*t)]; suT=su*T; % Basis signals
Esum= 0;
% 16-QAM signal waveforms corresponding to rectangular constellation
for i=1:L
   for l=1:L
      s(i,l,1)=2*i-L-1; s(i,l,2)=2*l-L-1; %In-phase/quadrature amplitude  
      Esum= Esum +s(i,l,1)^2 +s(i,l,2)^2;
      ss(L*(l-1)+i,:)=[ssc(i,:) sss(l,:)];
      sw(L*(l-1)+i,:)=s(i,l,1)*su(1,:)+s(i,l,2)*su(2,:); 
   end
end
Eav=Esum/M, Eav_s=2*(M-1)/3 % Eq.(7.5.4a): Average signal energy (A=1)
Es=2; % Energy of signal waveform
A=sqrt(Es/Eav); sw=A*sw; levels=A*[-(L-1):2:L-1];
SNRdBs=[1:15]; MaxIter=10000; % Range of SNRbdB and # of iterations 
for iter=1:length(SNRdBs)
   SNRbdB= SNRdBs(iter);  SNR=10^(SNRbdB/10);
   N0=2*(Es/b)/SNR;  sigma2=N0/2;  sgmsT=sqrt(sigma2/T);
   yr= zeros(2,LB);  nobe= 0; % Number of bit errors to be accumulated
   for k=1:MaxIter
      im= ceil(rand*L); in= ceil(rand*L);
      imn= (in-1)*L+im; % Index of signal to transmit
      s=ss(imn,:); % Data bits to transmit       
      for n=1:Ns   % Operation per symbol time
         wct= wcT*(n-1);  bp_noise= randn*cos(wct)-randn*sin(wct);
         rn= sw(imn,n) + sgmsT*bp_noise;
        yr= [yr(:,2:LB) suT(:,n)*rn]; % Multiplier
      end
      ycsk=sum(yr(:,LBN1:LB)'); % Sampled correlator output - DTR input
      %Detector(DTR)
      [dmin_i,mi]= min(abs(ycsk(1)-levels));
      [dmin_l,ml]= min(abs(ycsk(2)-levels));
      d= ss((ml-1)*L+mi,:); % Detected data bits
      nobe = nobe+sum(s~=d);  if nobe>100; break; end
   end
   pobe(iter)= nobe/(k*b);
end   
subplot(222), semilogy(SNRbdBt, pobet, 'k-', SNRdBs, pobe, 'b*')
title('Probability of Bit Error for 16-ary QAM Signaling')