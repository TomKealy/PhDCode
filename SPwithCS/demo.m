%% demo.m
% 
% This file contains a short demo that compares the DPSS version of CoSaMP
% to the DFT version of CoSaMP.
%
% Most recent change - 9/16/2011
%
% Copyright 2011, Mark Davenport, Michael Wakin
%
% This file is part of DPSS Approximation and Recovery Toolbox version 1.0.
%
%    DART is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    DART is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DART.  If not, see <http://www.gnu.org/licenses/>.

%% Initialize
clear all; close all;

addpath('../')

%% Set signal/algorithm parameters
N = 1024*4;
B = 1/256;
K = 5;

senseOpts.sampler = 'rdmdemod';  % Select random demodulator
senseOpts.pnseq = sign(randn(N,1)); % Choose random pn sequence
senseOpts.subsamp = 2^3; % Subsample by factor of 8

sigOpts.N = N; sigOpts.B = B; sigOpts.K = K;
sigOpts.type = 'randtones';  
sigOpts.TPB = 50; 

% Cosamp Opts
cosampOpts.dpssSize = N*B + 20;
cosampOpts.dpssWidth = N*B/2;
cosampOpts.N = N;
cosampOpts.maxiter = 25;
cosampOpts.tol = 1e-3;
cosampOpts.debias = 1;

%% Test Cosamp
[x supp] = generateSignal(sigOpts); % Generate signal
y = sense(x,0,senseOpts); % Compute measurements
[x_dpss ~] = dpssCosamp(y,senseOpts,K,B,cosampOpts); % Recover via Cosamp

%% Compare to fft-based approach
[x_fft supp_hat] = fftCosamp(y,senseOpts,K*25,cosampOpts);

%% Display results
disp(['DPSS Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x_dpss))) 'dB']);
disp(['FFT Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x_fft))) 'dB']);

figure(1)
plot([1:N]/N-.5,20*log10(abs(fftshift(fft(x)))),'Linewidth',2)
hold on
plot([1:N]/N-.5,20*log10(abs(fftshift(fft(x_dpss)))),'r','Linewidth',2)
plot([1:N]/N-.5,20*log10(abs(fftshift(fft(x_fft)))),'g','Linewidth',2)
set(gca,'LineWidth',2,'FontSize',13,'FontName','Times New Roman');
set(gca,'YLim',[0 100])

title('Frequency domain')
legend('Original','DPSS','FFT')