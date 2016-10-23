function plot_MOD(T,lt,msg,modul,demodul,How,detected,Bd,Ad)
% plots AM signals and their spectra
Fs=1/T; % Sampling Frequency/Period
t=[1:lt]*T; f =[-Fs/2: Fs/lt: Fs/2];  % Time/Freq. vector
M=fftshift(fft(msg)); 
M=[M M(1)]*T; % Spectrum of Message signal
Modul=fftshift(fft(modul)); 
Modul=[Modul Modul(1)]*T; % Spectrum of modulated signal
Y=fftshift(fft(demodul)); 
Y=[Y Y(1)]*T; % Spectrum of demodulated signal
subplot(421), plot(t,msg)
title('Message signal m(t)') 
subplot(422), plot(f,abs(M))
title('Spectrum of message')
subplot(423), plot(t,modul)
title([How ' modulated signal']) 
subplot(424), plot(f,abs(Modul))
title('Spectrum of modulated signal') 
subplot(425), plot(t,demodul)
title('Demodulated signal y(t)') 
subplot(426), plot(f,abs(Y))
title('Spectrum Y(f) of y(t)')
if nargin==9
  H=fftshift(fft(Bd,lt)./fft(Ad,lt));
  Hm=abs([H H(1)]); % Frequency Response of LPF
  hold on, plot(f,Hm,'r-')
end
if nargin>6
  Y_dtr=fftshift(fft(detected)); 
  Y_dtr=[Y_dtr Y_dtr(1)]*T; % Spectrum of detected signal
  subplot(427), plot(t,detected)
  title('Lowpass filtered output y_dtr(t)') 
  subplot(428), plot(f,abs(Y_dtr))
  title('Spectrum Y_dtr of y_dtr(t)')
end
