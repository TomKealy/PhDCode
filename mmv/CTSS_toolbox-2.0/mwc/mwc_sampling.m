function [y,Phi]=mwc_sampling(x,q,L,hb,ha,M)
% Modulated wideband converter sampling of a simulated sparse multiband signal.

% Usage: [y,Phi]=mwc_sampling(x,q,L,hb,ha,M)
% x: discrete input signal that approximates x(t)
% q: number of channels
% L: specifies rate and period of random +/- 1 sequence wrt Nyquist rate (1/Tp=W/L)
% hb, ha: FIR low-pass filter coefficients
% M: specifies sampling rate per channel (1/Ts=W/M)
% y: sampled output
% Phi: matrix of random +/- 1 (size L x q) 

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

Phi=pnmat(q,L); % generate random +/- 1 sequences
Phi_exten(:,:)=Phi(:,mod(0:length(x)-1,L)+1); % periodic extensions of Phi
y=zeros(q,length(x));
for i=1:q
 y(i,:)=x.*Phi_exten(i,:); % multiply input signal by random +/- 1 sequence 
end
y=filter(hb,ha,y.'); % filter demodulated signal with low pass filter
y=(downsample(y,M))'; % sample

function [Phi]=pnmat(M,m)
% Generates a M x m random matrix whose entries are drawn from {-1,+1} with equal probability.
% Usage: [Phi]=pnmat(M,m)

s=RandStream.getDefaultStream;
reset(s); % uncomment to reset random number generator and get repeatable results
Phi=sign(rand(s,[M,m])-0.5);  