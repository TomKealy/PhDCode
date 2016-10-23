%dc02f05.m
% generates white noise and computes its autocorrelation and DFT-PSD
clear, clf
rand('twister',5489); %return rand() to its default initial state
N=500;  % Length of Data Sequence
N_DFT=128;  % Size of DFT and length of Autocorrelation to be computed
phiw_av=zeros(1,N_DFT);  PHIw_av=zeros(1,N_DFT);
for iter=1:20	  % To take the ensemble average over 5 realizations
   w=rand(1,N)-1/2;%N i.i.d. r.v.s uniformly distributed over [-1/2,1/2]
   phiw= xcorr(w,w,'bias'); % autocorrelation of w[n]
   phiw= phiw(N-N_DFT/2+1:N+N_DFT/2);    
   PHIw= fftshift(fft(phiw)); % PSD of w[n]
   phiw_av= phiw_av+phiw;  % sum(ensemble average) of correlations
   PHIw_av= PHIw_av+PHIw; % sum(ensemble average) of PSD
   PHIwmag= abs([PHIw_av PHIw_av(1)]);
   subplot(311)
   plot([0:99], w(1:100)), title('w[n]')
   subplot(312)
   plot([-N_DFT/2+1:N_DFT/2], phiw_av/iter)
   title('Autocorrelation phix[n]')
   subplot(313)
   plot([-0.5: 1/N_DFT: 0.5], PHIwmag/iter)
   title('Power Spectral Density')  
   if iter<5, pause; end
end
