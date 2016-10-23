function [Xk,kk,xht,tt]=CTFS(x,P,N,ng)
if nargin<4, ng=221; end
if nargin<3, N=10; end
kk=[-N:N];  w0=2*pi/P; 
T=2*P; tt=[-T:T/100:T];
xejkw=[x '(t).*exp(-j*k*w0*t)'];
xejkwt=inline(xejkw,'t','k','w0');
tol=0.001;
for k=0:N
  X(k+1)=quadl(xejkwt,-P/2,P/2,tol,[],k,w0);
  %X(k+1)=quad(xejkwt,-P/2,P/2,tol,[],k,w0); %No problem with tol=0.0001
  %X(k+1)=smpsns(xejkwt,-P/2,P/2,100,k,w0);
end
Xk =[conj(X(N+1:-1:2)) X]; % to make the spectrum symmetric
X_mag= abs(Xk);  %Xph= angle(Xk);
k=1:N; jkw0t=j*k.'*w0*tt;
xht=(2*real(X(k+1)*exp(jkw0t))+X(1))/P;
xt=feval(x,tt); %original signal
if nargout<1
  subplot(ng), plot(tt,xt,'k-', tt,xht,'b:') 
  set(gca,'fontsize',9), %set(gca,'XLim',[-T T])
  axis([tt([1 end]) -0.2 1.2]), title('x(t) and xt=ICTFS(X(k)) up to Nth order')
  subplot(ng+1), stem(kk,X_mag,'MarkerSize',5,'LineWidth',1), set(gca,'fontsize',9)
  axis([kk([1 end]) -0.2 1.2]), title('CTFS Spectrum |X(k)|')
end
