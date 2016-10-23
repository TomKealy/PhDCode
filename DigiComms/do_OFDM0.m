%do_OFDM0.m
% simulates an OFDM system with Nfft=128-subcarrier DFT and 2^4=16-QAM
% where a frame (block) carries Nfft/2*4 (message) bits
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
Nfft=128; Nffth=Nfft/2;  % FFT size and its half
b=4; Nblock=Nffth*b; % 2^b QAM and Number of message bits per block
Ng=16; % # of samples of guard interval (consisting of prefix/suffix)
g=[0 1 0 0]; lc=length(g); % Channel response and its length
[gmax,d1]=max(g); d=d1-1; % Estimate of Channel delay
%d=d-1; % Wrong estimate of channel delay for STO effect
Nsym=Nfft+Ng; % Length of OFDM symbol
lx=Nsym*(1+ceil(lc/Nsym)); % Length of x
% This must be greater than or equal to the effective length of 
% channel response g
x=zeros(1,lx); % Buffer of transmitted sequence
ly=lx+Nsym; % Length of y >= x + one OFDM symbol size
y=zeros(1,ly); % Buffer of received vector
SNRbdB=10; SNRb=10^(SNRbdB/10); 
sigma=sqrt(1/(SNRb*Nblock))
prob_error(SNRbdB,'QAM',b,'BER')
nobe=0;  MaxIter=10000;  Target_no_of_error=100;
for iter=0:MaxIter
   msg=randint(1,Nblock); % Message vector
   S=QAM(msg,b);
   for i=2:Nffth,  X(i)=S(i);  X(Nfft-i+2)=X(i)';  end % Eq.(11.1.2)
   X(1) = real(S(1));  X(Nffth+1) = imag(S(1));
   xn = real(ifft(X,Nfft));  % Eq.(11.1.3)
   % To update the modulated sequence x to be transmitted   
   x = [x(1+Nsym:lx) xn(Nfft-Ng+1:Nfft) xn]; % Add CP   
   % To update the received sequence r
   y(1:ly-Nsym) = y(1+Nsym:ly);
   for n=Nsym:-1:1
      n1 = ly+1-n;
      y(n1) = sigma*randn;      
      for n2=1:min(lx+1-n,lc)
         y(n1) = y(n1) + x(lx+2-n-n2)*g(n2)/gmax; % Eq.(11.1.6)
      end
   end      
   y1=y(lx+d-Nfft+1:lx+d); % Removing/neglect the CP and
   %keps=0.1; y1=y1.*exp(j*2*pi*keps/Nfft*[0:Nfft-1]); % CFO effect
   Y = fft(y1,Nfft); % take the FFT
   Sd=Y; Sd(1)=real(Y(1))+j*real(Y(Nffth+1));
   subplot(222), hold on 
   plot(Sd+j*eps,'.')  % Received signal constellation
   % QAM demodulation
   detected = QAM_dem(Sd(1:Nffth),b);
   if iter>0
     nobe = nobe + sum(detected~=msg0);
     if nobe>Target_no_of_error, break; end 
   end
   msg0 = msg;
end
pobe = nobe/(iter*Nblock)