%do_uniform_quantization.m:  Uniform Quantization (Fig.4.1)
% gives boundary vector b, quantization level vector c, and MSQE
clear, clf
pdf='exp(-(x-m).^2/2/sigma^2)/sqrt(2*pi)/sigma'; %Gaussian pdf of x
xf=inline(['x.*' pdf],'x','m','sigma');
f=inline(pdf,'x','m','sigma');
m=0; sigma=1; % Mean and variance of the random variable x
b0=-3; bN=3; % Given least/greatest value of the random variable x
for N=5:6 % Number of quantization intervals
   delta=(bN-b0)/N;  b=b0+[0:N]*delta;
   msqe=0; % Mean-Square Quantization Error
   for i=1:N % Centroid of each interval
      tmp1=quad(xf,b(i),b(i+1),0.01,[],m,sigma);
      tmp2=quad(f,b(i),b(i+1),0.01,[],m,sigma);
      tmp=tmp1/tmp2;  c(i)=tmp;  % Eq.(4.1.1)
      x2f=inline(['(x-tmp).^2.*' pdf],'x','m','sigma','tmp');
      msqe=msqe+quad(x2f,b(i),b(i+1),0.01,[],m,sigma,tmp); % Eq.(4.1.2)
   end
   b,c % Resulting boundary vector and centroid vector
   x=b0+(bN-b0)/1000*[0:1000];
   y(find(x<b(1)))=c(1); % Left-most interval
   for i=1:N,  y(find(b(i)<=x&x<b(i+1)))=c(i);  end
   y(find(x>=b(N+1)))=c(N); % Right-most interval
   subplot(2,2,N-4), plot(x,y), hold on, grid on % Quantization graph
   m=0; sigma=1; fx=feval(f,x,m,sigma); plot(x,fx,'r:') ;
   msqe % Resulting MSQE
end
