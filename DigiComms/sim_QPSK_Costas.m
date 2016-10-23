%sim_QPSK_Costas.m
% simulates a QPSK system (in Fig.7.7) 
% with a carrier phase recovery using recursive Costas loop (Fig.8.7)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
KC=1; % 0/1 for without/with phase compensation
b=2; M=2^b; % Number of bits, Number of symbols (Modulation index)
Tb=1e-2; Ts=b*Tb; % Bit/Symbol time
Nb=40; Ns=b*Nb;   % Numbers of sample times in Tb and Ts
T=Ts/Ns; LB=4*Ns; LBN1=LB-Ns+1; % Sample time, Buffer size, ..
SNRbdBt=0:0.1:10; pobet=prob_error(SNRbdBt,'PSK',b,'bit'); % Eq.(7.3.8)
% QPSK signal waveforms
ss=[0 0; 0 1; 1 1; 1 0];  gs='>^<v';
wc=2*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T; wct=wc*t;
pi2=pi*2; phases=[0:M-1]*(pi2/M);
gam=0.01; % Stepsize of adjusting phase estimate
su=sqrt(2/Ts)*[cos(wct); -sin(wct)]; suT=su*T; %Basis signal Eq.(7.3.2)
thd=-pi/8; % Phase difference between XMTR and RCVR
Es=2;      % Energy of signal waveform representing the symbols, ..
for m=1:M  % PSK Signal waveforms
   sw(m,:)=sqrt(2*Es/Ts)*cos(wct+phases(m)+thd); % Eq.(7.3.1)
end
SNRbdBs=[1  5  9]; N_SNRbdBs=length(SNRbdBs); % Range of SNRdB
MaxIter=40000; % # of iterations for getting the error probability
Target_no_of_error=100;
for iter=1:N_SNRbdBs
   SNRbdB=SNRbdBs(iter);  SNRb=10^(SNRbdB/10);
   N0=2*(Es/b)/SNRb;  sigma2=N0/2;  sgmsT=sqrt(sigma2/T);
   yr=zeros(2,LB); % Multiplier output buffer
   thh(1)=0; % Phase offset estimate initialized
   nobe= 0; % Number of bit errors to be accumulated
   rand('twister',5489); randn('state',0); 
   for k=1:MaxIter
      i= ceil(rand*M); s=ss(i,:); % Signal index, Data bits to transmit
      for n=1:Ns % Operation per symbol time
         wct=(n-1)*wcT;
         rn=sw(i,n)+sgmsT*(randn*cos(wct)-randn*sin(wct));%Received signal
         yr= [yr(:,2:LB) suT(:,n)*rn]; % Multiplier output (buffer)
      end
      ycsk= sum(yr(:,LBN1:LB)'); % Correlator output samples - DTR input
      if iter==N_SNRbdBs&k<min(200,MaxIter) %Signal constellation diagram
        subplot(221), plot(ycsk(1),ycsk(2),gs(i)), hold on
      end
      ynejth=(ycsk(1)+j*ycsk(2)); % Not yet phase compensation
      if KC>0, ynejth=ynejth*exp(-j*thh(k)); end  % Phase compensation 
      th=angle(ynejth); % Compensated phase of received signal waveform
      if th<-pi/M, th=th+pi2; end
      [thmin,lmin]=min(abs(th-phases)); D=ss(lmin,:); % Detected bits
      thh(k+1)=thh(k)+gam*imag(exp(-j*phases(lmin))*ynejth); %Eq.(8.4.10)
      nobe=nobe+sum(s~=D);  if nobe>Target_no_of_error, break; end
   end
   pobe(iter)= nobe/(k*b); % BER
end   
subplot(222), semilogy(SNRbdBt,pobet,'k-', SNRbdBs,pobe,'b*')