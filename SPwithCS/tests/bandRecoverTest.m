%% Test bandRecover.m
% 
% Most recent change - 9/9/2011
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

if(0)
    senseOpts.sampler = 'rdmsample';  % Select random demodulator
    senseOpts.subsamp = 2^3; % Subsample by factor of 8
    rindex = randperm(N);
    senseOpts.idx = rindex(1:N/senseOpts.subsamp); % Select random indices
end

if(1)
    senseOpts.sampler = 'rdmdemod';  % Select random demodulator
    senseOpts.pnseq = sign(randn(N,1)); % Choose random pn sequence
    senseOpts.subsamp = 2^3; % Subsample by factor of 8
end

if(0)
    senseOpts.sampler = 'matrix';
    senseOpts.subsamp = 8;
    M = N/8;
    senseOpts.matrix = 1/sqrt(M)*randn(M,N);
end

sigOpts.N = N; sigOpts.B = B; sigOpts.K = K;
sigOpts.type = 'randtones';  
sigOpts.TPB = 50; 

%  DPSS KNOBS
dpssSize = N*B+20;
dpssWidth = N*B/2;    

%% Perform recovery
[x supp] = generateSignal(sigOpts); % Generate signal
modulator = generateModulator(N,B); % Generate modulator matrix
Psi = dpss(N,dpssWidth,dpssSize); % Construct base DPSS basis
y = sense(x,0,senseOpts); % Compute measurements
detectOpts.alg = 'blockOMP'; % Select detection algorithm
omega = bandDetect(x,2*K,Psi,Psi,modulator,detectOpts); % Estimate support of x

% Constructn incorrect support estimate 
KMax = floor(length(y)/size(Psi,2));
T = omega(1:min(KMax+1,end));
T = setdiff(T,supp(1));

% Recovery via correct support estimate
x_oracle = bandRecover(y,senseOpts,omega,Psi,modulator,norm(x)+1e-10);  

% Recovery via incorrect support estimate
x_hat = bandRecover(y,senseOpts,T,Psi,modulator,norm(x)+1e-10);  

%% Display results
disp(['N = ' num2str(N) ', M = ' num2str(N/senseOpts.subsamp) ', Dimension = ' num2str(K*dpssSize)])
disp(['Oracle SNR: ' num2str(20*log10(norm(x)/norm(x-x_oracle))) 'dB']);
disp(['Recovery SNR: ' num2str(20*log10(norm(x)/norm(x-x_hat))) 'dB']);

figure(1)
plot([1:N]/N-.5,abs(fftshift(fft(x))))
hold on
plot([1:N]/N-.5,abs(fftshift(fft(x_oracle))),'r')
title('Oracle recovery - Frequency domain')
legend('Original','Recovered')

figure(2)
plot([1:N]/N-.5,abs(fftshift(fft(x))))
hold on
plot([1:N]/N-.5,abs(fftshift(fft(x_hat))),'r')
title('Incorrect recovery - Frequency domain')
legend('Original','Recovered')