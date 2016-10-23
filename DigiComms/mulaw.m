function [y,xmax]=mulaw(x,mu,ymax)
if nargin<3, ymax=1; end
if nargin<2, mu=255; end
xmax=max(abs(x)); y=ymax*log(1+mu*abs(x/xmax))./log(1+mu).*sign(x);
