%sim_PCM.m
% simulates the PAM/PCM(Pulse Amplitude/Code Modulation) system
clear, clf
T=0.01; nmax=200; t=[0:nmax-1]*T;
b0=-1; bN=1; N=8; delta=(bN-b0)/N;  % Quantization interval
b=b0+[0:N]*delta; % Boundary vector
c=b(1:N)+delta/2; % Quantization level vector
Nb=ceil(log2(N)); N2=2^Nb;
if N2>N, c(N+1:N2)=c(N)*ones(1,N2.N); end
%code_table=fliplr(de2bi([0:N2-1]));
code_table=[0 1 0;0 1 1;0 0 1;0 0 0;1 0 0;1 0 1;1 1 1;1 1 0];%Gray code
noise_amp=0.51;
% PAM
signal_power=0;
for n=1:nmax
  x(n)=sin(6*(n-1)*T); 
  xq(n)=adc(x(n),b,c); %
  y(n)=xq(n)+noise_amp*(rand-0.5);
  yq(n)=adc(y(n),b,c); % Detection by RCVR
  signal_power=signal_power+xq(n)*xq(n);
end
tmp=ones(Nb,1)*xq;
transmitted_signal=tmp(:)';
Xts_PAM=abs(fftshift(fft(transmitted_signal)));
prob_of_symbol_error_with_PAM=sum(yq~=xq)/nmax
signal_power=signal_power/nmax % Power of PAM signal
%PCM
transmitted_signal=[]; % Serial bit pattern
signal_power=0; A=0.5;
for n=1:nmax
   [xq(n),code]=adc(x(n),b,c,code_table);
  for m=1:Nb % Nb=length(code)
      ym(m)=A*code(m)+noise_amp*(rand-0.5);
      ymq(m)=(ym(m)>A/2); % Detection by RCVR
  end    
   transmitted_signal=[transmitted_signal code(1:Nb)];
   y_dec(n)=c(strmatch(ymq,code_table));
  signal_power=signal_power+sum(code.^2);
end
Xts_PCM=abs(fftshift(fft(A*transmitted_signal)));
prob_of_symbol_error_with_PCM=sum(y_dec~=xq)/nmax
signal_power=A*A*signal_power/Nb/nmax % Power of PCM signal
subplot(311), plot(t,x,'b'), hold on, stairs(t,xq,'r')
w=2*pi/T/nmax/Nb*[-nmax*Nb/2:nmax*Nb/2]; % Frequency range
subplot(312), plot(w,[Xts_PAM Xts_PAM(1)],'b') % Spectrum of PAM signal
subplot(313), plot(w,[Xts_PCM Xts_PCM(1)],'b') % Spectrum of PCM signal
