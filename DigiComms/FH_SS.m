%FH_SS.m to simulate a digital communication system 
%  with coherent BFSK/fast FH_SS signaling in Fig.7.4.1 and Fig.10.7
clear, clf
m=4; % Number of stages in PN generator 
Nc=2^m-1; kh=3; % PN code length (Number of chips per message bit)
g=[1 4]; code1=png(g); % PN sequence
for i=1:Nc/kh
   wc(i)=code1(kh*(i-1)+1:kh*i)*2.^[kh-1:-1:0]'*2*Nc*pi;
end
b=1; M=2^b;
SNRbdBt=0:0.1:10;  SNRbt=10.^(SNRbdBt/10);
pobet=prob_error(SNRbdBt,'FSK',b,'BER');
Tb=1; Ts=b*Tb; Tc=Ts/Nc; % Bit/Symbol/Chip time 
N=20; T=Tc/N; % Sampling Period(subinterval in Tc)
NNc=N*Nc; N4=4*NNc; % Buffer size 
Es=2; sqEs=sqrt(Es); % Energy of signal waveform
ss=[0 1]; w=[8 10]*Nc*pi;
t=[0:3*N-1]*T; s2sT=sqrt(2/Ts);
SNRbdBs=[6 7]; PJ=0.01; MaxIter=1000;
for i_FH=1:2  % Frequency Hopping or not
   su0=zeros(Nc/3,3*N); sw=zeros(M,NNc);
   for m=1:M % PN-coded BPSK signal waveforms
      for nc=1:Nc/3
         if i_FH==1, su0(nc,:)=cos(w(m))*((nc-1)*3*N*T+t); % No FH
          else  su0(nc,:)=cos((w(m)+wc(nc))*((nc-1)*3*N*T+t)); % FH
         end  
         su(m,:)=s2sT*reshape(su0,1,NNc);
      end
   end  
   sw=sqEs*su; suT=su*T;
   for i_SNR=1:length(SNRbdBs) % Range of SNRdB
      SNRbdB=SNRbdBs(i_SNR); SNRb=10^(SNRbdB/10);
      N0=2*(Es/b)/SNRb; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
      Jamming_amp=sqrt(PJ/T); %sqrt(J0/Tc)
      r=zeros(1,N4);
      nobe=0; % Number of bit errors
      for k=1:MaxIter
         i=ceil(rand*M); s=ss(i); % Transmitted signal       
         for m=1:NNc   % Operation per symbol interval        
            tmp=w(i)*m*T; intfrnc=cos(tmp); % Interference
            r=[r(2:N4) sw(i,m)+sgmsT*randn+Jamming_amp*intfrnc]; 
            yre=suT*r(3*NNc+1:N4)'; % Correlator output
         end
         [yremax,lmax]=max(yre); D=ss(lmax); % Detector(DTR)
         nobe=nobe+sum(s~=D);
         if nobe>100, break; end
      end
      pobe(i_FH,i_SNR)=nobe/(k*b);
   end
end
pobe