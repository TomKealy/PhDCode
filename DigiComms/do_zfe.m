%do_zfe.m  to simulate ZFE (zero-forcing equalizer)
clear, clf
g=[0.027 0.0385 0.0588 0.1 0.2 0.5 1 0.5 0.2 0.1 0.0588 0.0385 0.027]; 
lg=length(g); %the channel response and its length(duration)
subplot(611), stem([0:lg-1],g)
title('Channel impulse response')
K=5; % Order (# of taps) of ZFE
for d=0:4 % to change the delay 
   c=zfe(g,K,d) % ZFE coefficients
   g_zfe=filter(c,1,g); % filtering the channel response with ZFE
   subplot(612+d), stem([0:lg-1],g_zfe) % plot the output of ZFE
end
