function [Xk,kk,xht,tt]=CTFT(x,D,ng)
if nargin<3, ng=221; end
N=100;
kk=[-N:N]; dw=10*pi/(2*N); 
T=5*D; tt=[-T:T/100:T];
xejkw=[x '(t).*exp(-j*k*w0*t)'];
xejkwt=inline(xejkw,'t','k','w0');
tol=0.001;
for k=0:N
  X(k+1)=quadl(xejkwt,-D,D,tol,[],k,dw);
  %X(k+1)=smpsns(xejkwt,-D,D,50,k,dw);
  if x(1:2)=='rD', Xw(k+1)=D*sinc(k*dw*D/2/pi); 
  elseif x(1:2)=='tr', Xw(k+1)=D*(sinc(k*dw*D/2/pi))^2;
  end
end
Xk =[conj(X(N+1:-1:2)) X]; % to make the spectrum symmetric
X_mag= abs(Xk);  %Xph= angle(Xk);
Xwk =[conj(Xw(N+1:-1:2)) Xw]; % to make the spectrum symmetric
Xw_mag= abs(Xwk);
%xt=zeros(size(tt));
%for k=1:N, xt=xt+X(k+1)*exp(j*k*dw*tt); end
%xt=(real(xt)+X(1)/2)*dw/pi;
k=1:N; jkwt=j*k.'*dw*tt;
xht=(real(X(k+1)*exp(jkwt))+X(1)/2)*dw/pi;
xt=feval(x,tt);
if nargout<1
  subplot(ng), plot(tt,xt,'k-', tt,xht,'b:')
  title('x(t) and xh(t)=ICTFT(X(w)) up to Nth order')
  axis([tt([1 end]) -0.5 1.5]), set(gca,'fontsize',9)
  subplot(ng+1), plot(kk,X_mag,'r-')
  if x(1:2)=='rD'|x(1:2)=='tr', hold on, plot(kk,Xw_mag,'k:'); end
  title('CTFT Spectrum |X(w)|'), set(gca,'fontsize',9)
end
