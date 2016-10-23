function [y,xmax]=Alaw(x,A,ymax)
if nargin<3, ymax=1; end
if nargin<2, A=87.56; end
A1=1/A;  xa=abs(x); xmax=max(xa); 
for n=1:length(x)
   if xa(n)<=A1, y(n)=ymax*A*x(n)/xmax/(1+log(A));
    else y(n)=ymax*(1+log(A*xa(n)/xmax))/(1+log(A))*sign(x(n));
  end
end 
