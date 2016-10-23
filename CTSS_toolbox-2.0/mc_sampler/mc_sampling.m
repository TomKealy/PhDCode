function [y,c]=mc_sampling(x,L,q,nyqover)
% This script implements the multi-coset sampling of a multiband signal.

% Usage: [y,c]=mc_sampling(x,L,q,nyqover)
% x: discrete input signal that approximates a continuous-time multiband signal
% L: parameter that determines the period of the nonuniform sampling (every 
%    L/W seconds q samples are collected); positive integer > q
% q: q samples collected every L/W seconds (q is also the number of channels)
% nyqover: controls the length of the simulated multiband signal x (conceptually, x is a vector 
%  formed from sampling a continuous-time multiband signal at a rate that is nyqover times more
%  than this Nyquist rate. nyqover is a positive integer.) 
% y: (sub)sampled output of multichannel system
% c: multi-coset sampling pattern

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

% construct sampling pattern (randomly chosen)
s=RandStream.getDefaultStream;
reset(s); % reset random number generator for repeatable results
pool=randperm(s,L);
c(:,1)=pool(1:q);

% sample according to multi-coset pattern  
y=zeros(floor(length(x)/(L*nyqover)),q);
for i=1:q
 y(:,i)=x(c(i)*nyqover:L*nyqover:floor(length(x)/(L*nyqover))*(L*nyqover)); % each column of y represents a time series of one channel 
end