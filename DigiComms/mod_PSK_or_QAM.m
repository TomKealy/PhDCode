function modulated=mod_PSK_or_QAM(in,b,Mod,Bin_or_Dec,Kg)
% Set Kg=1 to see the normal constellation diagram
if nargin<4&max(in)<2&length(in)==1, Bin_or_Dec='bin'; end
if max(in)>1, Bin_or_Dec='dec'; end
Lin= length(in); 
if lower(Bin_or_Dec(1))=='b' 
  % converts binary numbers into decimal numbers
  remi=rem(Lin,b); 
  if remi>0, in= [in zeros(1,b-tmp)]; Lin=Lin+1; end
  Lin = Lin/b; b1=b-1; 
  n=1; for m=1:Lin, in1(m)=bin2deci(in(n:n+b1)); n=n+b; end
  in = in1;
end  
M=2^b; 
if lower(Mod(1))=='p'
   constellation = exp(j*2*pi/M*[0:M-1]);
 else
   constellation = modulate(modem.qammod(M),[0:M-1]);
   constellation = constellation*modnorm(constellation,'avpow',1); % Normailized constellation.
end
dec_seq = gray2bin(in,Mod,M);
modulated = constellation(dec_seq+1);
% To see the constellation diagram
if nargout==0|(nargin>4&Kg>0)
  scatterplot(constellation) % Plot the normailized constellation.
  % Include text annotations that number the points.
  hold on; % Make sure the annotations go in the same figure.
  for k=1:length(constellation)
     text(real(constellation(k)),imag(constellation(k)),[' ' num2str(k-1)]);
  end
end  
%for n=1:lin
  %out(n)=exp(j*(phaseoffset+2*pi*gray2bin(in(n),b)/(b*2)));
  %out(n,:)=exp(j*pM2*gray2bin_my(in(n),b));  
%end
