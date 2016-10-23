function f=Rice_pdf(x,A,sigma)
%Rice density function 
sigma2= sigma^2; sgm2= 2*sigma2;
f=(x>=0).*x/sigma2.*??(x*A/sigma2).*exp(-(x.*x+A*A)/sgm2);
