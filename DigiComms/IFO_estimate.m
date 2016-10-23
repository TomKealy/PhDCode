function IFO_est=IFO_estimate(y,X,IFO_range)
% To estimate the IFO (Integral Frequency Offset)
% y: A received time-domain signal, supposedly containing ifft(PRS)
% X: (frequency-domain) Phase Reference Symbol
% IFO_range : Range of possible IFOs to be searched
M=3; Max=zeros(1,M); IFO_est=zeros(1,M); % 3 methods to estimate IFO
Nfft=length(X);  Y = fft(y,Nfft);
Ybar=Y.*conj(Y([Nfft 1:Nfft-1])); % Eq.(11.3.10) 
Xbar=X.*conj(X([Nfft 1:Nfft-1]));
for i=1:length(IFO_range)
   d=IFO_range(i);  YX = Y.*conj(rotate_r(X,d));
   Mag(1) = abs(sum(YX)); % Eq.(11.3.9a) APRS
   Mag(2) = max(abs(ifft(YX))); % Eq.(11.3.9b) ACIR 
   Mag(3) = abs(Ybar*rotate_r(Xbar,d)'); % Eq.(11.3.9c) AIDC
   for m=1:M
      if Mag(m)>=Max(m),  Max(m)=Mag(m); IFO_est(m)=d;  end
   end
end