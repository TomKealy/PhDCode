function z=corr_circular(x,y,N)
% Circular convolution z(n)=sum_m=0^N-1 x(m+n)*y(m)
if nargin<3, N=min(length(x),length(y)); end
x=x(1:N); y=y(1:N); x_circulated = x;
for n=1:N
   if n>1, x_circulated = [x_circulated(2:N) x_circulated(1)]; end
   z(n) = y*x_circulated'; % x_circulated*y';
end   
z= z/N;
