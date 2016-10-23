%
% Script for Generating Non-Negative Fourier Example
%
% solves the problem
% min ||x||_1 s.t. y = S*F*x, x >= 0 
%
% 
clear all;

% Signal
global p;
p = 2^9;
x = MakeBumps(p);x = x(:);

% Fourier transform and sampling
c = 8;
n = p/c;
global I;
I = 1:n;
z = fft(x);
y = z(I);

% Make System
b = [real(y);imag(y)];

% Solve
tic
[xhat, iters] = SolveIterSoftThresh('FastOp', b, p);           
toc

% Display
figure;
subplot(2,3,1); 
plot(x); 
axis([1 p 0 6]); 
title('signal');

subplot(2,3,2); 
plot(abs(fft(x))); 
axis([1 p 0 300]); 
title('Fourier transform of signal');

subplot(2,3,3); 
plot(abs(y)); 
axis([1 p 0 300]); 
title('Low frequencies');

subplot(2,3,4); 
plot(x,'--r','LineWidth',2); hold on;
ztag = z; ztag(setdiff(1:p,I)) = 0;
plot(real(ifft(ztag))); 
axis([1 p 0 6]); 
title('InvFourier of low freq')';

subplot(2,3,5); 
plot(x,'--r','LineWidth',2); hold on; 
plot(xhat); 
axis([1 p 0 6]); 
title('l1 reconstruction from low freq');

%
% Copyright (c) 2006. Iddo Drori
%

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
