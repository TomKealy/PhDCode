function [phi,lags]=xcorr_my(x,y,opt)
%computes the crosscorrelation of two vectors x and y by Eq.(2.1.61)
if nargin<3, opt='reg';  if nargin<2, y=x; end,  end
Nx=length(x); Ny=length(y); N=Nx+Ny-1; lags= [-(Ny-1):(Nx-1)];
for n=1:Nx-1
N1=min(Nx-n,Ny); m=1:N1; phi(n+Ny)=x(m+n)*y(m)'; 
if opt(1)=='u',  phi(n+Ny)= phi(n+Ny)/N1;  end
end
for n=1-Ny:0
N2=min(Nx,Nx+n); m=1:N2; phi(n+Ny)=x(m)*y(m-n)'; 
if opt(1)=='u',  phi(n+Ny)= phi(n+Ny)/N2;  end
end
if opt(1)=='b', phi= phi/N; 
 elseif opt(1)=='c', phi=phi/sqrt((x*x')*(y*y')); 
end
