function y= PSK_slicer(x,M)
% Converts received signal points into regular PSK points by its angle
y= zeros(size(x)); 
pi2=pi*2; %M=2^b; 
piM=pi/M; phs= pi2*[0:M-1]/M; 
MPSKs= exp(j*phs);
for n=1:length(x)
   ph= angle(x(n)); 
   if ph<-piM, ph= ph+pi2; end
   [emin,i]= min(abs(phs-ph));
   y(n)= MPSKs(i);
end
%plot(real(x),imag(x),'*')
