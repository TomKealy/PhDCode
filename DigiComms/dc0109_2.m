%dc0109_2.m
% Lowpass equivalent of Bandpass signal
clear, clf
ts =0.001; fs =1/ts; % Sampling Period/Frequency
N=512; N2=N/2;
t=[-(N2-1):N2]*ts; n1=[N2-40:N2+40]; t1=t(n1); % duration of signal
f1=200; w1=2*pi*f1; B=50; %center frequency & bandwidth of signal x(t)
f=[-fs/2:fs/N:fs/2]; % frequency vector
% Task 1
vt=3.5*sinc(2*B*t); pht=40*t;
x=vt.*cos(w1*t+angt); x1=x(n1); 
X=fftshift(fft(x)); % N=length(x)
subplot(421), plot(t1,x1) 
title('x(t)=3.5sinc(100*t)cos(2pi*200t+40*t)')
subplot(422), plot(f, abs([X X(1)]))
axis([-500 500 0 50]), title('Magnitude Spectrum of x(t)')
% Task 2
xa=hilbert(x); Xa=fftshift(fft(xa)); %xa(t)=x(t)+j x^(t)
subplot(423), plot(t1,imag(xa(n1)),'r'), hold on
th_tilt=0.2; % tilting angle of the real axis
costh=cos(th_tilt); sinth=sin(th_tilt);
xar=3*real(xa(n1)); % rough scaling
for n=1:length(t1)
xa_r(n)=t1(n)+xar(n)*costh/1000;  ya_r(n)=xar(n)*sinth;
end
plot(xa_r,ya_r), title('xa(t)=x(t)+j*xh(t)')
subplot(424), plot(f, abs([Xa Xa(1)])), 
title('Magnitude Spectrum of xa(t)')
% Task 3
xl=xa.*exp(-j*w1*t);  % Lowpass Equivalent Signal xl(t)
Xl=fftshift(fft(xl));
subplot(425), plot(t1,imag(xl(n1)),'r'), hold on
xlr=3*real(xl(n1)); % rough scaling
for n=1:length(t1)
   xl_r(n)=t1(n)+xlr(n)*costh/1000; yl_r(n)=xlr(n)*sinth;
end
plot(xl_r,yl_r), title('xl(t)=xa(t)*exp(-j*2pi*f1*t)')
subplot(426), plot(f, abs([Xl Xl(1)]))
title('Magnitude Spectrum of xl(t)')
% Task 4
env=abs(xa);  envl=abs(xl); % the envelope
subplot(427), plot(t1,env(n1), t1,envl(n1),'r', t1,vt(n1),'k:')
% Task 5
ph=principal_frequency(angle(xa)-w1*t); phl=angle(xl); % the angle
hold on, plot(t1,ph(n1), t1,phl(n1),'r', t1,pht(n1),'k:')
title('The envelope and angle')
Xh=fftshift(fft(imag(xa))); 
subplot(428), plot(f, abs([Xh Xh(1)]))
title('Magnitude Spectrum of hilbert transform x^(t)')
