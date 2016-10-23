%do_rcos1.m 
% Raised-cosine frequency response in Fig. P6.1
clear, clf
EPS=1e-8; pi2=pi*2; r=0.5; % Roll-off factor
Ts=1; T=Ts/4; % Symbol time and Sample time of filter output
Fd=1/Ts; Fs=Fd*(Ts/T); % Sampling frequencies of filter input/output
x=[1]; % An impulse input
% design a raised-cosine filter and implement it using rcosflt()
g2=rcosflt(x,Fd,Fs,'fir',r).'; N2=length(g2);
g3=rcosflt(x,Fd,Fs,'iir',r).'; N3=length(g3);
[gmax,im]=max(g2); 
N=im-1; N1=N*2+1; 
for n=0:N1-1
   if n~=N, t=(n-N)*T/Ts; else t=EPS; end;
   g1(n+1)=sinc(t).*cos(r*pi*t)./(1-(2*r*t).^2); % Eq.(6.1.7)
end
Nfft=256; Nfft2=Nfft/2;
G1=20*log10(abs(fftshift(fft(g1,Nfft)))); G1=[G1 G1(1)];
G2=20*log10(abs(fftshift(fft(g2,Nfft)))); G2=[G2 G2(1)];
G3=20*log10(abs(fftshift(fft(g3,Nfft)))); G3=[G3 G3(1)];
W0= 2*pi/Nfft/T; W=[-Nfft2:Nfft2]*W0;
subplot(221)
stem(0:N1-1,g1,'ko'), hold on
stem(0:N2-1,g2,'r.'), stem(0:N3-1,g3,'bx')
subplot(222), plot(W,G1, W,G2,'r', W,G3,'m'), hold on
% Theoretical raised-cosine frequency response 
r2=2*r; W1=(1-r)*pi/Ts; W2=(1+r)*pi/Ts; 
for k=1:length(W)
   absWk=abs(W(k));
   if absWk<W1, G0(k)=Ts*pi2;
    elseif absWk<W2
       G0(k)=Ts*(1+cos(Ts/r2*(absWk-W1)))*pi; % Eq.(6.1.3)
    else        G0(k)=1e-5;;
   end
end
G0= 20*log10(abs(G0));
plot(W,G0,'k') 
% design a raised-cosine filter 
%  and implement the specified fir by default
num = rcosine(Fd,Fs); g4 = rcosflt(x,Fd,Fs,'filter',num); 
% design a raised-cosine filter
%  and implement the specified iir/filter
delay=3; 
[num,den] = rcosine(Fd,Fs,'iir',r,delay); 
g5 = rcosflt(x,Fd,Fs,'iir/filter',num,den,delay);
[Max,im]=max(g5); N=im-1; N6=N*2+1; 
for n=0:N6-1
   if n~=N, t=(n-N)*T/Ts; else t=EPS; end;
   g6(n+1)=sinc(t).*cos(r*pi*t)./(1-(2*r*t).^2);
end
G4=20*log10(abs(fftshift(fft(g4.',Nfft)))); G4=[G4 G4(1)];
G5=20*log10(abs(fftshift(fft(g5.',Nfft)))); G5=[G5 G5(1)];
N4=length(g4); N5=length(g5); N0=(N5-N4)/2;
subplot(223) 
stem([0:N4-1]+N0, g4), hold on
stem([0:N5-1], g5), stem([0:N6-1], g6)
subplot(224), plot(W,G4,W,G5,'r',W,G0,'k')
