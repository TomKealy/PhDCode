%dc01e16.m 
% plot the CTFS and DFT spectra of a BFSK signal
clear, clf
Tb=0.001; f1=11; f2=13; D=Tb; P=2*D;
N=f1*f2; T=Tb/N; tt=[0:N-1]*T;
s1=cos(2*pi*f1/Tb*tt);  s2=cos(2*pi*f2/Tb*tt);
ttt=[tt Tb+tt];  ss=[s1 s2]; 
subplot(311)
plot(ttt,ss), axis([ttt([1 end]) -1.5 1.5])
title('A BFSK signal'), set(gca,'fontsize',9)
Nfft=length(ss); Nfft2=Nfft/2;
for k=-Nfft2:Nfft2
   Sk(k+Nfft2+1) = Sk_CTFS(k,D);
end
subplot(312)
ff=[-30:30]/(Nfft*T*1000);
stem(ff,abs(Sk(Nfft2-29:Nfft2+31))/T,'Markersize',5)
pause
for k=-Nfft2:Nfft2
   SSk(k+Nfft2+1)=Sk_CtFS(k,D)+Sk_CtFS(k-Nfft,D)+Sk_CtFS(k+Nfft,D)...
      +Sk_CtFS(k-2*Nfft,D)+Sk_CtFS(k+2*Nfft,D)...
      +Sk_CtFS(k-3*Nfft,D)+Sk_CtFS(k+3*Nfft,D);
end
stem(ff,abs(SSk(Nfft2-29:Nfft2+31))/T,'s','Markersize',5)
title('DFT: S(k)=1/T*\Sigma S_{k+mN} (The scaled sum of infinite shifted CTFS)')
xlabel('frequency[kHz]'), set(gca,'fontsize',9,'xtick',[-13 -11 0 11 13])
hold on, pause
S = fftshift(fft(ss));  % fftshift(fft(ss(N-Nfft2+1:N+Nfft2))); 
Smag=abs([S S(1)]); % DFT spectrum
stem(ff,Smag(Nfft2-29:Nfft2+31),'rx','Markersize',5)
pause, legend('Scaled sum of CTFSs','2*11*13=286point-DFT',0)
%title('Comparison of the CtFS and DFT spectra')
%set(gca,'fontsize',9)
err = norm(Smag(Nfft2-29:Nfft2+31)-abs(SSk(Nfft2-29:Nfft2+31))/T)