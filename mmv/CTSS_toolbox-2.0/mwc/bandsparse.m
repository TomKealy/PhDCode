function [x,centerfreq,t]=bandsparse(W,B,K,L,type,T)
% Script generates approximations to simple (and idealized) continuous-time wideband 
% signals that have K occupied frequency bands.
                                                               
% NOTE: In some cases, parameter values are set manually instead of random choices.

% Usage: [x,centerfreq,t]=bandsparse(W,B,K,L,type,T)
% W: Nyquist frequency (bandlimit of the multiband signal= W/2)
% B: bandwidth of each occupied frequency band (Hz)
% K: number of occupied bands (positive integer)
% L: parameter that determines the period of the nonuniform sampling (every L/W seconds q samples 
%    are collected); positive integer > q
% type: specifies the type of multiband signal (options are 'modsinc' (modulated sinc waveforms), 'filterwhitenoise' (white
%       noise filtered through a notch filter), and 'chirp' (multiband chirp signal)
% T: observation interval [0,T) sec
% x: simulated multiband signal
% centerfreq: center frequencies of occupied bands
% t: vector of discrete time instances over which band sparse signal is defined

% Copyright (c) 2010 by Michael A Lexa, Mike E Davies, John S Thompson
% Institute of Digital Communications, University of Edinburgh
%
% This work is licenced under the Creative Commons Attribution 3.0
% Unported License. To view a copy of this licence, visit
% http://creativecommons.org/licenses/by/3.0/ or send a letter to
% Creative Commons, 171 Second Street, Suite 300, San Francisco,
% California 94105, USA.

switch type

case 'modsinc'
Tend=floor(ceil(T/(1/W))/L)*L/W;
t=(0:1/W:Tend-1/W); % time instances over which simulated signal is defined

%centerfreq=rand(1,K)*W-W/2; % randomly generate center frequencies of occupied bands
centerfreq=[-150 20 350]; % use to manually set centerfreq; comment if randomly generating delays
%delays=rand(1,K)*T; % uncomment to randomly generate time delays of sinc pulses
delays=[2.5 8.0 15.0]; % use to manually set time offsets of sinc pulses; comment if randomly generating delays
tt=(-Tend:1/W:Tend-1/W);
x=zeros(1,length(t)); % initialize multiband signal vector

for n=1:K % generate K time-domain sinc pulses and modulate each of them to their respective center frequencies
 sincpul=sqrt(B)*sinc(B*(tt-delays(n))); % time domain sinc pulses
 x=x+sincpul(floor(Tend*W)+1:end).*exp(1i*2*pi*centerfreq(n)*t); % modulate pulses 
end


case 'filterwhitenoise' % occupied bands are filtered white Gaussian noise 
Tend=floor(ceil(T/(1/W))/L)*L/W;
t=(0:1/W:Tend-1/W); % time instances over which simulated signal is defined

%centerfreq=rand(1,K)*W/2; % uncomment to randomly generate center frequencies of occupied bands
centerfreq=[20 200 320]; % use to manually set centerfreq; comment if randomly generating center frequencies
x=zeros(1,length(t)); % initialize multiband signal vector
for n=1:K
 f=[(centerfreq(n)-(B/4))/(W/2) (centerfreq(n)+(B/4))/(W/2)];
 hb=fir1(200,f,'bandpass'); % window based FIR design
 %s=RandStream.getDefaultStream; reset(s); % uncomment for repeatable results
 x=x+filter(hb,1,randn(1,length(x))); % filter white Gaussian noise
end


case 'chirp' % creates linear chirp within each occupied band
Tend=floor(ceil(T/(1/W))/L)*L/W;
t=(0:1/W:Tend-1/W); % time instances over which simulated signal is defined

%centerfreq=rand(1,K)*W/2; % randomly generate center frequencies of occupied bands
centerfreq=[250 350 100]; % uncomment to manually set centerfreq 
x=zeros(1,length(t)); % initialize multiband signal vector
delays=[2 8 15]; % set chirp offsets manually
durations=[2 1 2]; % set chirp durations manually
for n=1:K
 chirpsig=zeros(1,length(t));
 chirppul=chirp(t(1:durations(n)*W),centerfreq(n),t(durations(n)*W),centerfreq(n)+B); % generate real chirp pulse
 chirpsig(delays(n)*W:(delays(n)+durations(n))*W-1)=chirppul;
 x=x+chirpsig;
end

end % switch
