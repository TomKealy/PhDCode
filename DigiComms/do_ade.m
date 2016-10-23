%do_ade.m  to simulate Adaptive equalizer
clear, clf
g=[0.027 0.0385 0.0588 0.1 0.2 0.5 1 0.5 0.2 0.1 0.0588 0.0385 0.027]; 
lg=length(g); % Channel response and its length(duration)
Na=500; % Length of training period
K=5;    % Order of MMSE EQLZR
[gmax,n01]=max(g); n0=n01-1; % channel delay 
d=n0+floor(K/2); % Expected delay of channel & equalizer (6.2.12)
c_ADE= zeros(1,K); % initialize the adaptive equalizer coefficients
ch_input= zeros(1,lg); % initialize the channel input buffer
sigma=0.2;  % Amplitude of noise
delta=0.05; % Step size
for m=1:Na % Training period
   a(m)= (rand>0.5)*2-1; % Training sequence to be available to RCVR
   ch_input= [a(m) ch_input(1:lg-1)]; % channel input
   y(m)= ch_input*g' + sigma*(rand-0.5);
   if m>max(K-1,d)
     [c_ADE,se(m)]=ade(c_ADE,y(m:-1:m-K+1),a(m-d),delta);
   end
end
gn = g+sigma*(rand(size(g))-0.5);
g_ADE =filter(c_ADE,1,gn); % ADE output to (channel response+noise)
subplot(414), stem([0:lg-1],g_ADE), title('ADE response')
