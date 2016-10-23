function [x,X,freq]=tonesparse(t,Tx,N,K)
% Generates a spectrally sparse multitone signal. 

% Usage: [x,X,freq]=tonesparse(t,Tx,N,K)
% t: time interval (vector) over which signal is observed (sec)
% Tx: length of observation interval 
% N: length of simulated input signal
% K: the number of nonzero FS coefficients/frequencies comprising x
% x: simulated spectrally-sparse time-domain multitone signal
% X: K-sparse vector of FS coefficients; nonzero elements are randomly chosen
% freq: support vector (multiples of fundamental frequency) 

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

n=-N/2:N/2-1; % range of possible frequency indices (multiples of fundamental frequency)
s=RandStream.getDefaultStream;reset(s);
freq=(n(ceil(rand(s,[K,1])*length(n))))'; % randomly pick K indices (determines support of multitone signal)  
%s=RandStream.getDefaultStream; reset(s); % uncomment to reset random number generator for repeatable results
%XX=rand(s,[1,length(freq)])*100; % uncomment to randomly generate amplitudes of FS coefficients (magnitudes in [0,100])	
XX=10*ones(1,1); % set all amplitudes to 10
x=(1/N)*sum(diag(XX)*exp(1i*(2*pi)/Tx*freq*t),1)'; % create multitone signal
X=zeros(N,1);
X(floor(freq)+N/2+1)=XX; % assign nonzero frequency amplitudes