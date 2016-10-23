%sim_biorthogonal.m 
% simulates an M=2^b-ary communication system with bi-orthogonal signal waveforms in Fig.5.18
clear, clf
b=2; M=2^b;  % # of bits/symbols
Tb=1; Nb=20; T=Tb/Nb; % Bit interval time Tb and Sample time T
Ts=b*Tb; Ns=b*Nb; % Symbol time Ts and Number of samples Ns per symbol time
LB=4*Ns; LB1=LB-Ns+1; % Length of biffer
Es=2; sqEs=sqrt(Es); 
su=[? ? ? ?;? ? -? -?]; su=[su; -su];
su=ones(Ns/M,1)*reshape(su',1,M*M);
su=reshape(reshape(su,1,Ns*M)',Ns,M)'/sqrt(Ts); suT=su*T;
sw=sqEs*su; % Bi-orthogonal signal waveforms
ss=[0 0;0 ?;? 0;? ?]; % Set of M-ary symbola represented in binary numbers
SNRbdBs=[1:10]; MaxIter=20000; % Range of SNRbdB and Number of iterations
for iter=1:length(SNRbdBs)
   SNRb=10^(SNRbdBs(iter)/10); N0=2*(Es/b)/SNRb;  sgmsT=sqrt(N0/2/T);
   st= zeros(1,LB);  y= zeros(M,LB);  r= zeros(1,LB);  yr= zeros(M,LB);
   nobe= 0; % Initialize the number of bit errors to be accumulated
   for k=1:MaxIter
      i= ceil(rand*M); s=ss(i,:); % Transmitted signal             
      for m=1:Ns % Operation per symbol interval
         r = [r(2:LB) sw(i,m)+sgmsT*randn]; % Received signal
         yre = suT(1:M/2,:)*r(LB1:LB)'; % Sampled correlator output
      end
      [ymax,imax]=max(abs(???)); % Detector(DTR)
      if yre(imax)>=0, d=ss(imax,:);  else d= ss(imax+???,:);  end
      nobe = nobe+sum(d~=s);  if nobe>100, break; end
   end
   pobe(iter) = nobe/(k*b); % Bit error probability
end
SNRbdBt=0:0.1:10; SNRbt=10.^(SNRbdBt/10); 
tol=1e-10; tmp=M/2/(M-1)/sqrt(2*pi);
f5260_s='(1-(1-2*Q(x)).^(2^(b-1)-1)).*exp(-(x-sqrt(??????)).^2/2)'; %Eq.(5.2.60)
f5260 = inline(f5260_s,'x','SNRb','b');
for i=1:length(SNRbt), pobet(i)=tmp*quadl(f5260,0,40,tol,[],SNRbt(i),b); end
semilogy(SNRbdBt, pobet, 'k-', SNRbdBs, pobe, 'b*')
title('BER for M-ary Bi-Orthogonal Signaling')