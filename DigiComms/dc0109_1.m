%dc0109_1.m
clear, clf
N=64; Ts=1/16; % sampling period
t= [-N/2+0.5:N/2-0.5]*Ts; % time vector without t=0
h= 1/pi./t; % the impulse response of the Hilbert transformer
Nh= length(t); Nh2= floor(Nh/2);
h_= fliplr(h); % h(-t): the time-reversed version of h(t)
wc1= pi; wc2= 2*pi; % the frequency of an input signal
Nfft=N; Nfft2= Nfft/2; Nbuf=2*Nfft;
ww=[-Nfft2:Nfft2]*(2*pi/Nfft); % frequency vector
tt= zeros(1,Nbuf);  x_buf= zeros(1,Nbuf);  xh_buf= zeros(1,Nbuf);
for n=0:Nbuf-1
   tn=n*Ts; tt=[tt(2:end) tn]; 
   x_buf = [x_buf(2:end) sin(wc2*tn)]; %sin(wc*t)]; 
   xh_buf = [xh_buf(2:end) h_*x_buf(end-Nh+1:end).'*Ts];
end
subplot(431)
plot(tt,x_buf), title('x(t)')
set(gca,'Xlim',tt([1 end]),'Ylim',[-1.2 1.2], 'fontsize',9)
X= fftshift(fft(x_buf,Nfft)); X(find(abs(X)<1e-10))=0;
X=[X X(1)]; % X(w): spectrum of x(t)
subplot(432), plot(ww,abs(X)), title('|X(w)|')
subplot(433), plot(ww,angle(X)), title('<X(w)')
subplot(434), plot(t(1:Nh2),h(1:Nh2), t(Nh2+1:end),h(Nh2+1:end))
% circular-shift for causality
h_(1:N/2)=h(N/2+1:end); h_(Nfft-N/2+1:Nfft)=h(1:N/2); 
H= fftshift(fft(h_,Nfft)); H=[H H(1)];
H1= -j*sign(ww)/Ts;
subplot(435), plot(ww,abs(H), ww,abs(H1),'r:')
subplot(436), plot(ww,angle(H)-0.5*ww, ww,angle(H1),'r:') 
% To advance the delayed response of the causal Hilbert transformer
xh_advanced = xh_buf(Nh2+1:end); 
Xh_1 = -j*sign(ww).*X;
xh_1 = real(ifft(fftshift(Xh_1(1:end-1)),Nfft));
xa = hilbert(x_buf); % analytic signal using hilbert()
xh_2 = imag(xa); % Hilbert transform using hilbert()
discrepancy_t = norm(xh_1-xh_2(1:length(xh_1)))/norm(xh_2)
subplot(437) 
plot(tt(1:end-Nh2),xh_advanced, tt(1:length(xh_1)),xh_1,'m.') 
hold on, plot(tt,xh_2,'r:'), title('xh(t)')
Xh= fftshift(fft(xh_advanced,Nfft)); Xh= [Xh Xh(1)]; % Xh(w)
subplot(438), plot(ww,abs(Xh), ww,abs(Xh_1),'r:') 
subplot(439), plot(ww,angle(Xh_1),'r:'), title('<Xh(w)')
subplot(4,3,10)
th=0.8; % the tilting angle of the real axis
for n=1:length(tt)
plot(tt(n)+[0 real(xa(n))],[0 imag(xa(n))+real(xa(n))*tan(th)])
hold on
end   
xa=hilbert(x_buf); % xa(t): the analytic signal
Xa=fftshift(fft(xa,Nfft)); Xa=[Xa Xa(1)]; Xa1=2*(ww>0).*X;
discrepancy_w=norm(Xa-Xa1)
subplot(4,3,11), plot(ww,abs(Xa), ww,abs(Xa1),'r:') 
subplot(4,3,12), plot(ww,angle(Xa1),'r:'), title('<Xa(W)')
