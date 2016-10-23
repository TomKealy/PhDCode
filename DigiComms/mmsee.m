function  c=mmsee(y,a,K,d)
%Input:  y= Output of channel to a (training) data sequence a
%        a= Message data, K=# of equalizer taps
%        d= Expected channel&equalizer delay 
%Output: c= Equalizer coefficients
ly=length(y); phi_y=xcorr(y,'unbiased'); phi_ay=xcorr(a,y,'unbiased');   
for i=1:K
   Ry(i,:)=phi_y(ly-i+1:ly-i+K); % Eq.(6.2.10)
   rya(i)=phi_ay(ly+i-1-d); % Eq.(6.2.11)
end
c=rya/Ry; % Eq.(6.2.9)
