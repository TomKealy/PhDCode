%FH_SS2.m
% simulates a digital communication system 
%  with coherent BFSK-slow FH_SS signaling in Fig.10.7.2(b)
clear, clf
b=1; M=2^b; kh=3;
alpha=1/(2^kh); % depending on the jammer's strategy
Tb=1; Ts=b*Tb; Ncs=2; Tc=Ncs*Ts; % Bit/Symbol/Chip time 
N=48; T=Ts/N; % Sampling Period (subinterval in Tc)
NNcs=N*Ncs; % Buffer size 
for i=1:2^kh,  wc(i)=(i-1)*M*2*pi;  end
Es=2; sqEs=sqrt(Es); % Energy of signal waveform
ss=[0 1]; w([1:M])=2*[1:M]*pi;  g_symbol='xo';
t=[0:N-1]*T;  s2sT=sqrt(2/Ts);
SNRbdBs=[6]; MaxIter=3000; % # of iterations 
for i_FH=1:2
   sw=zeros(M,NNcs);
   for m=1:M % PN-coded BPSK signal waveforms
      for nc=1:2^kh
         if i_FH==1, su(m+(nc-1)*M,:)=s2sT*cos(w(m)*t); % No FH
          else  su(m+(nc-1)*M,:)=s2sT*cos((w(m)+wc(nc))*t); % FH
         end  
      end
   end  
   sw=sqEs*su;  suT=su*T;
   for i_SNR=1:length(SNRbdBs) % Range of SNRdB
      SNRbdB=SNRbdBs(i_SNR); SNRb=10^(SNRbdB/10); J0=2*(Es/b)/SNRb; sigma=sqrt(J0/Tb*(2^kh)); %Pj=J0*W=J0*(2/Tb)*2^kh
      r=zeros(1,N);
      nobe=0; % Number of bit errors
      for k=1:MaxIter     
         ic=floor(rand*2^kh);
         for k1=1:Ncs
            i= ceil(rand*M); s=ss(i); % Transmitted signal
            for m=0:N-1   % Operation per symbol interval
               mT=m*T; tmp=(w(1)+w(2))/2*mT; intfrnc=sin(tmp);
               r=[r(2:N) sw(i+M*ic,m+1)+sigma*intfrnc]; % Received signal
               yre=suT(([1:M]+M*ic)',:)*r'; % Correlator output 
            end
            [yremax,lmax]=max(yre);  D=ss(lmax); % Detector(DTR)
            nobe=nobe+sum(s~=D);  if nobe>100, break; end
         end
      end  
      pobe(i_FH,i_SNR)=nobe/(k*b);
   end       
   semilogy(SNRbdBs,pobe(i_FH,:),g_symbol(i_FH))
end
pobe