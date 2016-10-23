%sim_QAM_carrier_recovery.m 
% simulates a QAM system with Decision-feedback Carrier Recovery in Fig. 8.9.
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use 
clear, clf
KC=2; % 0/1/2 for no phase offset/no carrier recovery/carrier recovery
b=4; M=2^b; L=2^(b/2);
SNRbdBt=0:0.1:15; pobet=prob_error(SNRbdBt,'QAM',b,'bit'); %Eq.(7.5.6)
Tb=0.025; Ts=b*Tb; % Bit/Symbol time 
Nb=20; Ns=b*Nb; % Numbers of sample times in Tb and Ts
T=Ts/Ns; LB=Ns*4; LBNs=LB-Ns; % Sample time, Buffer size
% 16-QAM signal waveforms
ssc=[0 0; 0 1; 1 1; 1 0]; sss=ssc; % In_phase/Quadrature data bits
wc=12*pi/Ts; t=[0:Ns-1]*T; s2sT=sqrt(2/Ts); % Carrier frequency, time
su= s2sT*[cos(wc*t); -sin(wc*t)]; % Basis signal waveforms
for m=1:L
  for n=1:L
     s(m,n,1)=2*m-L-1; s(m,n,2)=2*n-L-1;   
     ss(L*(n-1)+m,:)=[ssc(m,:) sss(n,:)]; % Data bits per symbol
     sw(L*(n-1)+m,:)=s(m,n,1)*su(1,:)+s(m,n,2)*su(2,:); 
  end
end
Es=2; % Energy of signal waveform
Eav=2*(M-1)/3; A=sqrt(Es/Eav); sw=A*sw; levels=A*[-(L-1):2:L-1].';
zeta=0.707; wn=pi/2000/Ts; % Damping ratio and Natural frequency of LF
% Parameters of Loop Filter to be designed
T2= max(2*zeta/wn,T/pi), T1=10*T2, KvT= 4*wn^2*T1*T;  
% Discrete-time LF through BLT with prewarping
[bLF,aLF]=bilinear([T2 1],[T1 0],1/T,1/T2/2/pi)
% Still, no general method of tuning the LF and VCO is clear to us. 
SNRbdBs=[3 8 13]; 
MaxIter=4e4; % Number of iterations for getting the error probability
Target_no_of_error=100; 
for iter=1:length(SNRbdBs)  
   SNRbdB=SNRbdBs(iter);  SNRb=10^(SNRbdB/10);
   N0=2*(Es/b)/SNRb; sigma2=N0/2; sigmaT=sqrt(sigma2/T);
   r= zeros(1,LB); yr= zeros(2,LB); ycs= zeros(2,LB);
   ve=zeros(1,LB); vLF=zeros(1,LB); 
   thk=0; Ak=1; nobe= 0; % Number of bit errors to be accumulated
   rand('twister',5489); randn('state',0); 
   for k=1:MaxIter
      im= ceil(rand*L); in= ceil(rand*L); imn= (in-1)*L+im;
      s=ss(imn,:); % Data bits to transmit
      if k<5, Nd=0; thh=0;  
       elseif KC>0, Nd=1;  
      end
      thd=-wc*Nd*T; % True phase offset
      for n=1:Ns % Operation per symbol time
         wct= wc*(n-1)*T;
         bp_noise= randn*cos(wct)-randn*sin(wct);
         r=[r(2:LB) sw(imn,n)+sigmaT*bp_noise]; 
         yr=[yr(:,2:LB) s2sT*[cos(wct+thh); -sin(wct+thh)]*r(LB-Nd)];  
         ycs=[ycs(:,2:LB) ycs(:,LB)+yr(:,LB)-yr(:,LBNs)];
         ve=[ve(2:LB) [-sin(thk) cos(thk)]/Ak*yr(:,LBNs)]; 
         vLF=[vLF(2:LB) vLF(LB)+bLF*ve(LB:-1:LB-1)']; %Phase shift by thk
         thh= thh +KvT*vLF(LB); % Phase offset estimate
         if KC<2, thh=0; end % No carrier recovery
      end
      ycsk=ycs(:,LB)*T; 
      yck=ycsk(1); ysk=ycsk(2);
      if iter==length(SNRbdBs)&k<200 % Signal constellation diagram
        subplot(341+KC), hold on, plot(ycsk(1),ycsk(2),'.')
      end
      %Detector(DTR)
      [dmin,lmn]=min(abs([ycsk(1)-levels ycsk(2)-levels])); 
      lm=lmn(1); ln=lmn(2); % Detected (in-phase/quad) signal indices
      thk=atan2(levels(ln),levels(lm));
      D= ss((ln-1)*L+lm,:); % Detected data bits
      nobe=nobe+sum(s~=D); % Number of bit errors
      if nobe>Target_no_of_error, break; end
   end
   pobe(iter)= nobe/(k*b);
end   
subplot(345+KC)
semilogy(SNRbdBt,pobet,'k-', SNRbdBs,pobe,'b*')
title('BER for (16-ary) QAM Signaling')
