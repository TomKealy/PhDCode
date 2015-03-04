%Author: Osama Ullah Khan,
%        Phd Student, University of Michigan-Ann Arbor.
%        Email: oukhan@umich.edu
%        Version: 1.0
%
%This code demonstrate compressive sensing example. In this
%example the signal is sparse in frequency domain and random samples
%are taken in time domain.

close all;
clear all;

%setup path for the subdirectories of l1magic
path(path, 'C:\MATLAB7\l1magic-1.1\Optimization');
path(path, 'C:\MATLAB7\l1magic-1.1\Data');


%length of the signal
N=1024;

%Number of random observations to take
K=256;

%Discrete frequency of two sinusoids in the input signal
k1=29;
k2=100;

n=0:N-1;

%Sparse signal in frequency domain.
x=sin(2*pi*(k1/N)*n)+sin(2*pi*(k2/N)*n);

figure;
subplot(2,1,1);
plot(x)
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal,1024 samples with two different frequency sinsuoids');

xf=fft(x);

xfmag=10*log10(abs(xf));

subplot(2,1,2);
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Frequency domain, 1024 coefficients with 4-non zero coefficients');

%creating dft matrix
B=dftmtx(N);
Binv=inv(B);

%Taking DFT of the signal
xf=B*x';

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=Binv(q(1:K),:);

%taking random time measurements
y=(A*xf);

%Calculating Initial guess
x0=A'*y;

%Running the recovery Algorithm
% tic
% xp=l1eq_pd(x0,A,[],y,1e-5);
% toc

cvx_begin
    variable x(N) complex;
    minimise (norm(xf,1))
    subject to
        A*xf == y
cvx_end

%recovered signal in time domain
xprec=real(Binv*xf);

figure;
subplot(2,1,1)
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal, Discrete Fourier Transform');

subplot(2,1,2)
plot(abs(xf),'r')
grid on;
xlabel('Samples');
ylabel('Amplitude');
title(sprintf('Recovered Signal, Discrete Fourier Transform sampled with %d samples',K));

figure;
subplot(2,1,1);
plot(xf)
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal,1024 samples with two different frequency sinsuoids');

subplot(2,1,2)
plot(xprec,'r')
grid on;
xlabel('Samples');
ylabel('Amplitude');
title(sprintf('Recovered Signal in Time Domain'));