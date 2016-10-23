%do_dfe.m  to simulate Decision Feedback Equalizer(DFE)
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
clear, clf
g=[0.027 0.0385 0.0588 0.1 0.2 0.5 1 0.5 0.2 0.1 0.0588 0.0385 0.027]; 
lg=length(g); % Channel response and its length(duration)
[gmax,n01]=max(g);  
n0=n01-1; % Dragging time or peak time of channel response 
K=5; % the order of EQLZR
Lt=500; Ls=5000; % Lengths of training and simulation periods 
sigma=0.01; % Amplitude of noise
delta=0.005; % Step size
d=n0+floor(K/2); % Expected delay of channel and equalizer (6.2.12)
c=zeros(1,K);  b=zeros(1,K); % initialize the DFE coefficients
% initialize the history buffers for channel input/output, DTR output
ch_input= zeros(1,lg); ys= zeros(1,K); ds= zeros(1,K);
neb=0; % initialize the number of error bits
for n=1:Ls
   a(n)= (rand>0.5)*2-1; % Message data
   ch_input= [a(n) ch_input(1:lg-1)]; % Channel input history
   y= ch_input*g'+sigma*(rand-0.5); % Corrupted channel output
   ys= [y  ys(1:K-1)]; % Channel output history
   q(n)=c*ys'-b*ds'; % DTR input by Eq.(6.2.16)
   D=(q(n)>0)*2-1;  % DTR output
   e=0; % initialize the discrepancy between DTR input and output
   if d<n&n<=Lt
     e=a(n-d)-q(n); % Eq.(6.2.18) for training period
     elseif n>Lt % after training period
       e=D-q(n); % Discrepancy between DTR input and output
       neb=neb+(D~=a(n-d)); % Number of error bits
   end  
   [c,b]=dfe(c,b,ys,ds,e,delta); % DFE
   ds= [D ds(1:K-1)]; % DTR output history
   se(n)=e^2; % Squared discrepancy between DTR input and output
end  
pbe=neb/(Ls-Lt) % Error probability
c,b
plot(se)