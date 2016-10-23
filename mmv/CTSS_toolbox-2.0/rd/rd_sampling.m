function [y,S]=rd_sampling(x,L,h,N,M)
% This script randomly demodulates the input signal x and uniformly samples the result.

% Usage: [y,S]=rd_sampling(x,L,h,N,M)
% x: input multitone signal(approximates continuous-time multitone signal)
% L: specifies period of random +/- 1 sequence wrt Nyquist rate (1/Tp=W/L)
% h: impulse response of ideal integrator
% N: length of simulated input signal
% M: the ratio M/N specifies rate of sampling process relative to Nyquist rate 
% y: sampled output
% S: vector of random +/- 1's 

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

S=pnmat(L,1); % generate random +/- 1 sequences
y=x.*S; % multiply input signal by random +/- 1 sequence 
y=conv(h,y); % filter demodulated signal with ideal integrator
y=downsample(y(N/M:N),N/M); % sample


function [Phi]=pnmat(M,m)
% Generates a M x m random matrix whose entries are drawn from {-1,+1} with equal probability.
% Usage: [Phi]=pnmat(M,m)

Phi=sign(rand(M,m)-0.5);     