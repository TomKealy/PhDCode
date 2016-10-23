%do_S_QDPSK.m 
% simulates the pi/4-shifted QDPSK signaling (Table P7.5)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
b=2; M=2^b; M2=M*2;
SNRbdBt=0:0.1:10;   SNRbt=10.^(SNRbdBt/10);
pobet=prob_error(SNRbdBt,'DPSK',b,'bit');
Tb=1; Ts=b*Tb; % Bit/symbol time 
Nb=16; Ns=b*Nb; % Numbers of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time and Buffer size
Es=2; % Energy of signal waveform
% QDPSK signal waveforms
ss=[0 0; 0 1; 1 1; 1 0]; 
wc=8*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T;  nd=1;
for m=1:M2,  sw(m,:)=sqrt(2*Es/Ts)*cos(wc*t+(m-1)*pi/M);  end
su= sqrt(2/Ts)*[cos(wc*t); -sin(wc*t)]; suT=su*T;
SNRdBs =[1:3:10]; MaxIter=10000; %Range of SNRbdB, # of iterations
for iter=1:length(SNRdBs) 
   SNRbdB=SNRdBs(iter);  SNRb=10^(SNRbdB/10);
   N0=2*(Es/b)/SNRb; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
   sws=zeros(1,LB); yr=zeros(2,LB);
   nobe=0; % Number of bit errors to be accumulated
   is0=1; % Initial signal index
   th0=1; % Initial guess (possibly wrong)
   for k=1:MaxIter
      i= ceil(rand*M); s=ss(i,:); % Data bits to transmit
      is= mod(is0+?????,M2)+1; % Signal to transmit (Table P7.5)
      for n=1:Ns % Operation per symbol time
        sws=[sws(2:LB) sw(is,n)];
         wct=wcT*(n-1); bp_noise= randn*cos(wct)-randn*sin(wct);
         rn= sws(end-nd) + sgmsT*bp_noise;
         yr=[yr(:,2:LB) suT(:,n)*rn]; % Multiplier
      end
      ycsk=sum(yr(:,LBN1:LB)')'; % Sampled correlator output
      %Detector(DTR)
      th=atan2(ycsk(2),ycsk(1)); dth=th-th0;
      if dth<0, dth=dth+2*pi; end
      [themin,lmin]=min(abs(dth-[?????]*2*pi/M));
      d= ss(lmin,:); % Detected data bits
      nobe= nobe+sum(s~=d);  if nobe>100, break; end
      is0=is; th0=th; % update the previous signal and theta
   end
   pobe(iter)= nobe/(k*b);
end   
semilogy(SNRbdBt,pobet,'k', SNRdBs,pobe,'b*')
