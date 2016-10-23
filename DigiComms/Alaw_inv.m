function x=Alaw_inv(y,A,xmax)
if nargin<3, xmax=1; end
if nargin<2, A=87.56; end
ya=abs(y); yA=max(ya)/(1+log(A));
for n=1:length(y)
   if ya(n)<=yA, x(n)=xmax*ya(n)/yA/A*sign(y(n));
    else x(n)=xmax*exp(ya(n)/yA-1)/A*sign(y(n));
  end
end 
