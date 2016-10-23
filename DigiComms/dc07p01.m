%dc07p01.m 
% See if a linear convolution of Gaussian noises 
% with a unit signal waveform is another Gaussian noise 
clear, clf
K=10000; % # of iterations for getting the error probability
Ts=1; N=40; T=Ts/N; % Symbol time and Sample time
N4=N*4; % Buffer size of correlator
wc=10*pi/Ts;  t=[0:N-1]*T;  wct=wc*t;
su=sqrt(2/Ts)*cos(wct); % Unit energy signal 
signal_power=su*su'*T % Signal energy
sigma2=2; sigma=sqrt(sigma2); sqT=sqrt(T);
noise= zeros(1,N4); % Noise buffer
for k=1:K
   for n=1:N % Operation per symbol time
      noise0= sigma*randn; 
      noise=[noise(2:N4) noise0/sqT]; % Bandpass noise
      noise1= su*noise(3*N+1:N4)'*T;
   end
   noise0s(k)=noise0; noise1s(k)=noise1;
end 
phi0=xcorr1(noise0s);  phi1=xcorr1(noise1s);
plot(phi0), hold on, pause, plot(phi1,'r')
