function x=mulaw_inv(y,mu,xmax)
if nargin<3, xmax=1; ymax=1;
 else ymax=max(abs(y));
end
if nargin<2, mu=255; end
x=xmax.*(((1+mu).^(abs(y)/ymax)-1)/mu).*sign(y);
