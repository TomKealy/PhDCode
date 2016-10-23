%DS_SS.m
% To simulate a BPSK/DS-SS communication system
clear, clf 
g=[2 5]; code1=png(g); % A PN sequence
g=[2 3 4 5]; code2=png(g); % Another PN sequence
Nc=length(code1); % PN code length (Number of chips per message bit)
b=1; M=2^b; 
Tb=1; Ts=b*Tb; Tc=Tb/Nc; % Bit/Symbol/Chip time 
N=8; T=Tc/N; % Sampling Period (subinterval in Tc)
NNc=N*Nc;  Ns=b*NNc; % Number of samples in one symbol
Es=2; sqEs=sqrt(Es); % Energy of signal waveform
ss=[0 1]; % Set of symbols 
wc=4*Nc*pi/Tb;  t=[0:Ns-1]*T;  s2sT=sqrt(2/Ts);
SNRbdBt=0:0.1:10; pobet= prob_error(SNRbdBt,'PSK',b,'BER');
SNRbdBs=[6 7]; PJ=0.01;  MaxIter=20000; 
for i_DS=1:2 % DS coding or not
   if i_DS==1, code1s=ones(1,Ns); % No DS 
    else  code1s=reshape(ones(N,1)*(2*code1-1),1,Ns); % PN-DS coding
   end      
   for m=1:M % (Possibly PN-coded) BPSK signal waveforms
     sw(m,:)=sqEs*s2sT*(ss(m)*2-1)*code1s.*cos(wc*t); % 0/1 -> -1/1
   end
   suT= s2sT*T*[code1s.*cos(wc*t); -code1s.*sin(wc*t)];  
   for i_SNR=1:length(SNRbdBs) % Range of SNRdB
      SNRbdB=SNRbdBs(i_SNR);  SNRb=10^(SNRbdB/10);
      sigma2=(Es/b)/SNRb;  sgmsT=sqrt(sigma2/T);
      Jamming_amp=sqrt(PJ/T);
      y= zeros(1,Ns); % Buffer size amounting to 1 symbol time
      nobe= 0; % Number of bit errors
      for k=1:MaxIter
         i= ceil(rand*M); s=ss(i); %transmitted signal       
         for m1=1:b*Nc % Operation per symbol interval Ts
            for m=1:N % Operation per chip interval Tc
               mi=(m1-1)*N+m; tmp=wc*((k-1)*Ns+mi-1)*T; 
               intfrnc=cos(tmp); % Interference
               rNs= sw(i,mi) +sgmsT*randn +Jamming_amp*intfrnc;
               y(1,:)=[y(1,2:Ns) suT(1,mi)*rNs]; % Inside correlator  
            end
         end
         ycsk=sum(y'); % Sampled correlator output - DTR input
         lmin=(ycsk(1)>0)+1; D=ss(lmin); % Detector(DTR)
         nobe=nobe+sum(s~=D); if nobe>100, break; end
      end   
      pobe(i_DS,i_SNR)=nobe/(k*b);
   end
end
fprintf('\n SNRbdB= %5.2f   %5.2f (Jamming power=%5.2f)', SNRbdBs,PJ)
fprintf('\n    BER= %6.4f  %6.4f (Without PN coding)', pobe(1,:))
fprintf('\n    BER= %6.4f  %6.4f (With PN coding)\n', pobe(2,:))
semilogy(SNRbdBt,pobet, SNRbdBs,pobe(1,:),'x', SNRbdBs,pobe(2,:),'o')