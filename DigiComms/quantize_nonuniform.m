%quantize_nonuniform.m (Fig.4.2)
% Input : Least/Greatest value of x, the # of quantization intervals
% Output: boundary vector b, quantization level vector c,
%             mean-square quantization error(MSQE)
clear, clf
b0=-3; bN=3; % Given least/greatest of the random variable x
N=6; % the number of quantization intervals
%Initialize the quantization intervals to uniform ones
delta=(bN-b0)/N;  b=b0+[0:N]*delta
pdf='exp(-(x-m).^2/2/sigma^2)/sqrt(2*pi)/sigma'; %Gaussian pdf of x
xf=inline(['x.*' pdf],'x','m','sigma');
f=inline(pdf,'x','m','sigma');
m=0; sigma=1; % Mean and variance of the random variable x
msqe=1000;  MaxIter=40;
for iter=1:MaxIter
  msqe0=msqe; 
  msqe=0; %Mean-Square Quantization Error
  for i=1:N %centroid of each interval
     tmp1=quad(xf,b(i),b(i+1),0.01,[],m,sigma);
     tmp2=quad(f,b(i),b(i+1),0.01,[],m,sigma);
     c(i)=tmp1/tmp2; 
     x2f=inline(['(x-ci).^2.*' pdf],'x','m','sigma','ci');
     msqe=msqe+quad(x2f,b(i),b(i+1),0.01,[],m,sigma,c(i));
  end
  for i=2:N %new boundaries between intervals
    b(i)=(c(i-1)+c(i))/2;
  end
  subplot(211), 
  plot(b,iter*ones(size(b)),'bx'), hold on
  if msqe>=0.99*msqe0, break;  end
end
b,c
x=b0+[0:1000]*(bN-b0)/1000;
for i=1:N
   y(find(b(i)<=x&x<=b(i+1)))=c(i);
end
subplot(212), plot(x,y)
y_or_n=input('Do you want to compare with another result?','s');
if y_or_n=='y'
  b(2)=-2.16; b(3)=-1.51; b(4)=-.98; b(5)=-.48;  
  b(10)=2.16; b(9)=1.51; b(8)=.98; b(7)=.48;  
  msqe1=0;
  for i=1:N %centroid of each interval
    tmp1=quad(xf,b(i),b(i+1),0.01,[],0,1);
    tmp2=quad(f,b(i),b(i+1),0.01,[],0,1);
    c(i)=tmp1/tmp2;  
    x2f=inline(['(x-ci).^2.*' pdf],'x','m','sigma','ci');
    msqe1=msqe1+quadl(x2f,b(i),b(i+1),0.01,[],0,1,c(i));
  end
  msqe, msqe1
  subplot(211), plot(b,(iter+1)*ones(size(b)),'rx')
end