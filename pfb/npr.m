% Near perfect reconstruction polyphase filter bank demo
%
% (c) 2007-2010 Wessel Lubberhuizen

clear all;
close all;

% number of taps per channel 
L = 128;
% number of channels
N = 32;
K=11.4;

display('designing prototype filter...');
c=npr_coeff(N,L,K);

display('generating a test signal...');
M = 128; % number of slices

% generate a linear chirp as a test signal.
% matlab's own chirp function has too much phase noise
% so we use our own version!
t=(0:M*N-1)/(M*N);
dphi=t;
phi=zeros(size(dphi));
for i=2:length(dphi);
    phi(i) = mod(phi(i-1)+dphi(i-1),1);
end
x = exp(-sqrt(-1)*2*pi*phi);
length(x)
% add some white noise if you like
%x = awgn(x,200);

positions = randi(L,[1,3]);%generate random spikes for signal

edges =  [0, 500, 1500, 2500, 3500, 4096];
levels = [0, 0, 1, 0, 0, 1];
idxs = zeros(1, M)  ;
idxs(edges(1: end-1)+1) = 1 ;
g = levels(cumsum(idxs)+1);

g(length(g):M*N)=0;

Tx = x;

% run it through the npr filterbank
display('processing...');
time1 = cputime;
y=npr_analysis(c,Tx);
z=npr_synthesis(c,y);
time2=cputime;
display(sprintf('processing rate = %f kSamples / second',length(x)/(time2-time1)/1E3));

% compare the input and the output
% delay=N*(L-1)/2;
% padding = zeros(1,delay);
% z = z(1+delay:length(z));
% x = x(1:length(x)-delay);

display(sprintf('average reconstruction error = %f dB',20*log10(norm(z-Tx,2)/norm(Tx,2))));

figure();
subplot(2,1,1);
plot((Tx));
title('input / output signals');
xlabel('sample');
ylabel('signal value');
grid on;
subplot(2,1,2);
plot(abs(z'));
xlabel('time (sample)');
ylabel('error value');
grid on;

figure();
plot(20*log10(abs([Tx' z'-Tx'])));
title('reconstruction error');
xlabel('time (sample)');
ylabel('power (dB)');
grid on;

figure();

img=20*log10(abs(y(1:size(y,1),1:size(y,2))));
imagesc(img);
colorbar;
title('spectrogram (dB)');
xlabel('time (slice)');
ylabel('frequency (channel)');
