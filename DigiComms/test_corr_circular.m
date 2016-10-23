%test_corr_circular.m
clear, clf
Nt=160; T=1/32; t=[0:Nt-1]*T; P=1; N=P/T;
X = sin(2*pi/P*t);  y = cos(2*pi/P*t);
z = corr_????????(x,y,N); % Circular correlation of x(t) and y(t)
subplot(311), plot(t,x), title('x(t)=sin(2*pi*t)')
set(gca,'fontsize',9, 'YLim',[-1.1  1.1])
subplot(312), plot(t,y), title('y(t)=cos(2*pi*t)')
set(gca,'fontsize',9, 'YLim',[-1.1  1.1])
subplot(313), plot(t,repmat(z,1,Nt/N))
title('z(t)=circular correlation of x(t) and y(t)')
set(gca,'fontsize',9, 'YLim',[-1.1  1.1])
