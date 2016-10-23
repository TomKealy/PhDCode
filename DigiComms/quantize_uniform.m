%quantize_uniform.m (Fig.4.1)
% gives boundary vector b, quantization level vector c,
% mean-square quantization error(MSQE)
clear, clf
%Gaussian probability density function of x
pdf='exp(-(x-m).^2/2/sigma^2)/sqrt(2*pi)/sigma';
%pdf='exp(-(x-m).^2/2/sigma^2)';
xf=inline(['x.*' pdf],'x','m','sigma');
f=inline(pdf,'x','m','sigma');
m=0; sigma=1; % Mean and variance of the random variable x
b0=-3; bN=3; % Given least/greatest value of the random variable x
for N=5:6 % the number of quantization intervals
  delta=(bN-b0)/N;  b=b0+[0:N]*delta;
  msqe=0; %Mean-Square Quantization Error
  for i=1:N %centroid of each interval
     tmp1=quad(xf,b(i),b(i+1),0.01,[],m,sigma);
     tmp2=quad(f,b(i),b(i+1),0.01,[],m,sigma);
     tmp=tmp1/tmp2;  c(i)=tmp;
     x2f=inline(['(x-tmp).^2.*' pdf],'x','m','sigma','tmp');
     msqe=msqe+quad(x2f,b(i),b(i+1),0.01,[],m,sigma,tmp);
  end
  b,c
  x=b0+[0:1000]*(bN-b0)/1000; N1=N+1;
  %ind0=find(x<b(1)); x(ind0)=b(1)*ones(size(ind0)); %left-most interval
  %indN=find(x>b(N1)); x(indN)=b(N1)*ones(size(indN)); %right-most interval
  y(find(x<b(1)))=c(1); y(find(x>=b(N1)))=c(N);
  for i=1:N
     y(find(b(i)<=x&x<b(i+1)))=c(i);
  end
  subplot(2,2,N-4), plot(x,y) %quantization graph
  hold on, grid on
  fx=feval(f,x,m,sigma); %probability density ftn
  plot(x,fx,'r:')
  axis([-3 3 -3 3])
  msqe
end