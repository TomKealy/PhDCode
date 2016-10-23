function demodulated=dem_PSK_or_QAM(in,b,Mod,Bin_or_Dec,Kg)
Lin= length(in); demodulated = [];
M=2^b; 
if lower(Mod(1))=='p', constellation = exp(j*2*pi/M*[0:M-1]);
 else
  constellation = modulate(modem.qammod(M),[0:M-1]);
  constellation = constellation*modnorm(constellation,'avpow',1);
end
for n=1:Lin, [minv,ii(n)] = min(abs(in(n)-constellation)); end
demodulated = bin2gray(ii-1,Mod,M);
if nargin<4|lower(Bin_or_Dec(1))=='b', demodulated = deci2bin(demodulated); end
if nargout==0|(nargin>4&Kg>0)
  scatterplot(constellation) % Plot the received signal constellation.
  hold on, plot(real(in),imag(in),'*', constellation,'ro')
end

