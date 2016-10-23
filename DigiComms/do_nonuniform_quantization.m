%do_nonuniform_quantization.m
% Lloyd-Max Nonuniform Quantization (Fig.4.2)
clear, clf
b0=-3; bN=3; % Given least/greatest of the random variable x
N=6; % Number of quantization intervals
% Initialize the quantization intervals to uniform ones
delta=(bN-b0)/N;  b=b0+[0:N]*delta
pdf='exp(-(x-m).^2/2/sigma^2)/sqrt(2*pi)/sigma'; %Gaussian pdf of x
xf=inline(['x.*' pdf],'x','m','sigma');
f=inline(pdf,'x','m','sigma');
m=0; sigma=1; % Mean and variance of the random variable x
msqe=1000;  MaxIter=40;
for iter=1:MaxIter
  msqe0=msqe; % Previous MSQE
  msqe=0; % Mean-Square Quantization Error
  for i=1:N 
     tmp1=quad(xf,b(i),b(i+1),0.01,[],m,sigma);
     tmp2=quad(f,b(i),b(i+1),0.01,[],m,sigma);
     c(i)=tmp1/tmp2;  % Centroid of each interval by (4.1.5)
     x2f=inline(['(x-ci).^2.*' pdf],'x','m','sigma','ci');
     msqe=msqe+quadl(x2f,b(i),b(i+1),0.01,[],m,sigma,c(i)); %(4.1.2)
  end
  for i=2:N % New boundaries between intervals
    b(i)=(c(i-1)+c(i))/2; %(4.1.6)
  end
  subplot(211), 
  plot(b,iter*ones(size(b)),'bx'), hold on
  if msqe>=0.99*msqe0, break;  end % Check for convergence
end
b,c
x=b0+[0:1000]*(bN-b0)/1000;
for i=1:N,  y(find(b(i)<=x&x<=b(i+1)))=c(i);  end
subplot(212), plot(x,y)
