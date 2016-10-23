function  c=zfe(g,K,d)
%Input:  g= Channel response, i.e., Impulse response of channel 
%        K= Order (# of taps) of ZFE
%        d= Expected delay of equalized response
%Output: c= Equalizer tap coefficients
if nargin<3, d=ceil(K/2); end 
[gmax,n01]=max(g);  n0=n01-1;
g(end+1:n0+K)=0;
for m=1:K, G(m,:)=g(m+n0:-1:m+n01-K);  end
delta=zeros(K,1); delta(d+1)=1; 
c=G\delta; % Eq.(6.2.5)
