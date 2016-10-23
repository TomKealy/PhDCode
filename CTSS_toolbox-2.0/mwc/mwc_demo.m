%% Sampling sparse mulitband signals with a multi-channel random demodulator
% This top-level script (m-file) demonstrates the sampling and recovery of sparse
% multiband signals using the Modulated Wideband Converter (MWC). The MWC was originally 
% proposed by M. Mishali and Y. Eldar in "From Theory to Practice: Sub-Nyquist Sampling
% of Sparse Wideband Analog Signals", Selected Topics in Signal Processing, IEEE Journal
% of, vol. 4, no. 2, pp. 375-391, Apr. 2010.

% When the script is executed MATLAB will generate a discrete multiband signal 
% (simulating a continuous one), sample it, and recover the original signal from 
% its samples. The modeling and system parameters are set in this script. For a 
% complete description of the parameters, refer to the technical report,
% "Sampling Sparse Multiband Signals with a Modulated Wideband Converter", accompanying 
% the CTSS Sampling Toolbox or the above reference.

% Note that this script does not implement a digital compensation filter that has been
% shown to counteract the effects of non-ideal low pass filters. Note also that the
% sub-routine mwc_recovery.m uses the m-file RunOMP_Unnormalized.m written by M. Mishali. 
% It is available at http://www.technion.ac.il/~moshiko/software.html.

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

clear

%% Parameters for input signal and modulated wideband converter
W=1e3; % Nyquist frequency (Hz) (signal bandwidth= W/2 Hz)
B=20; % maximum bandwidth of each occupied frequency band (Hz)
K=3; % number of occupied frequency bands present the input signal x
T=1; % duration of observation interval in sec
q=10; % number of channels
M=50; % specifies sampling rate per channel (1/Ts=W/M)
L=50; % L specifies period of p(t) wrt the Nyquist period (Tp=L/W); setting L=M means period of p(t) equals sampling period
hfilterorder=200; % set FIR filter order (should be even)
hb=fir1(hfilterorder,1/(2*L)); ha=1; % coefficients characterizing window-based FIR LPF (cut-off freq=W/2M Hz)
SNR=100; % signal-to-noise ratio (in dB)

%% Parameter checks
if M>L 
 error('M must be less than or equal to L to prevent destructive aliasing');
end
if q>M
 error('q must be less than or equal to M to ensure sub-Nyquist sampling');
end

%% Simulated signal generation, (sub)sampling, support recovery, signal reconstruction
[x,centerfreq,t]=bandsparse(W,B,K,L,'filterwhitenoise',T); % generate simulated sparse multiband signal (Options: modsinc, filterwhitenoise, chirp)
input_sigP=norm(x)^2/length(x); % compute signal power
input_noiseP=(1/SNR)*input_sigP; %(1/(10^(SNR/10)))*input_sigP; % determine required noise power (variance) wrt input signal power
s=RandStream.getDefaultStream;
reset(s); % uncomment to reset random number generator for repeatable results
w=sqrt(input_noiseP)*(randn(s,size(x))+1i*randn(s,size(x))); % generate white Gaussian noise with required variance 
xw=x;%+w; % add noise

[yw,Phi]=mwc_sampling(xw,q,L,hb,ha,M); % MWC sampling
[centerfreqw_hat,xw_hat,Rw,Psi,A]=mwc_recovery(yw,L,W,Phi,M,t);  % recover support and reconstruct signal

%% Results
fprintf('Sampling window duration %3.3f sec \n',T);
fprintf('Nyquist rate %3.3f Hz\n',W);
fprintf('System sampling rate %3.3f Hz\n',(q*W)/M);
fprintf('Subsampling factor %3.3f\n',q/M);
fprintf('Input SNR %3.3f dB\n',10*log10(input_sigP./input_noiseP));
fprintf('True center frequencies of bands with bandwidth %3.2f Hz\n',B);
fprintf('%3.2f\n',sort(centerfreq)); % display true center frequencies of active bands
fprintf('Center frequencies of declared active bands (spectral resolution %3.2f Hz)\n',W/L);
fprintf('%3.2f\n',sort(centerfreqw_hat)); % display recovered center frequencies 
fprintf('Average squared error %2.10f\n',sum((real(xw(1:end-hfilterorder/2))-real(xw_hat(hfilterorder/2+1:end))).^2)*(1/(length(xw)-1))); % display average squared error
fprintf('Column rank of A matrix: %2i\n',rank(A)); % display rank of A
fprintf('Column rank of R matrix: %2i\n',rank(Rw)); % display rank of R

figure % plot fft of original signal x and reconstructed signal x_hat
f= (-length(x)/2:length(x)/2-1)./length(x) * W;
h0=subplot(311); plot(f,abs(fftshift(fft(x))),'Marker','none'); title('Fourier Spectrum (original signal)'); 
xlabel('Hz'); ylabel('Magnitude')
h1=subplot(312); plot(f,abs(fftshift(fft(xw))),'Marker','none'); title('Fourier Spectrum (noisy signal)'); 
xlabel('Hz'); ylabel('Magnitude')
h2=subplot(313); plot(f,abs(fftshift(fft(xw_hat))),'Marker','none'); title('Fourier Spectrum (reconstructed signal)');
xlabel('Hz'); ylabel('Magnitude')
linkaxes([h0 h1 h2])

figure % plot original multiband signal x and reconstructed signal x_hat
h3=subplot(411); plot(t,real(x),'Marker','none');
xlabel('Seconds'); title('Original signal (real part)');
h4=subplot(412); plot(t,real(xw),'Marker','none');
xlabel('Seconds'); title('Noisy signal (real part)');
h5=subplot(413); plot(t(1:end-(hfilterorder/2+1)),real(xw_hat(hfilterorder/2+1:end-1)),'r','Marker','none'); % shift in time related to delay in FIR low pass filter used in sampling
xlabel('Seconds'); title('Reconstructed signal (real part)');
h6=subplot(414); stem(t(1:end-(hfilterorder/2+1)),real(xw(1:end-(hfilterorder/2+1)))-real(xw_hat(hfilterorder/2+1:end-1)),'Marker','none'); 
xlabel('Seconds'); title('Difference Signal (noisy-reconstructed)');
linkaxes([h3 h4 h5 h6])

figure % graphically display key matrices
subplot(221); imagesc(Phi); axis equal; axis tight; title('\Phi matrix');
subplot(222); imagesc(abs(Psi)); axis equal; axis tight; title('\Psi matrix (magnitude)');
subplot(223); imagesc(abs(A)); axis equal; axis tight; title('A matrix (magnitude)');
subplot(224); imagesc(abs(Rw)); axis equal; axis tight; title('R matrix (magnitude)');