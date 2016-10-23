%sim_QPSK_Costas_earlylate.m 
% to simulate a QPSK system (Fig.7.7) 
% with a carrier phase recovery scheme using the recursive Costas loop
% and a timing recovery scheme using the recursive Early-Late algorithm
%Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
KC=1; 
% KC=0 for simulation with none of phase compensation & timing recovery
% KC=1 for simulation with phase compensation & no timing recovery
% KC=2 for simulation with phase compensation & timing recovery
gamma=0.01; % Stepsize of adjusting the phase estimate
eth=0.5; % Error threshold
b=2; M=2^b; % Number of bits and symbols
Tb=1e-5; Ts=b*Tb; % Bit/Symbol time
Nb=40; Ns=b*Nb; % Number of sample times in Tb and Ts
T=Ts/Ns;        % Sampling interval
LB=4*Ns; Ns_2=Ns/2; % Buffer size  
g=[0.1 -0.2 0.6 0.9 1 0.9 0.6 -0.2 0.1]; g=g/norm(g)^2;
Ng=length(g); zi_XMTR=zeros(1,Ng-1); 
ndd=1; % Detection delay due to XMTR and channel
% QPSK signal waveforms
ss=[0 0; 0 1; 1 1; 1 0]; 
wc=2*pi/Ts; wcT=wc*T; t=[0:Ns-1]*T; wct=wc*t;
pi2=pi*2;  phases=[0:M-1]*(pi2/M);
su=sqrt(2/Ts)*[cos(wct);-sin(wct)]; suT=su*T; % Basis signal Eq.(7.3.2)
thd=-pi/10; % Phase offset between XMTR and RCVR
Es=2; % Energy of signal waveform representing the symbols
for m=1:M  % PSK signal waveforms
   sw(m,:)=sqrt(2*Es/Ts)*cos(wct+phases(m)+thd); % PSK signal (7.3.1)
end
interval=2*Ns+[1:Ns-1]; % Regular integration interval [k-2,k-1]*Ts
interval1=interval+Ns_2; % Later integration interval [k-3/2,k-1/2]*Ts
SNRbdB=9; SNRb=10^(SNRbdB/10); N0=2*(Es/b)/SNRb; sgmsT=sqrt(N0/2/T);
pobet= prob_error(SNRbdB,'PSK',b,'bit') % Theoretical BER
yr= zeros(2,LB); % Multiplier output buffer
s0s=zeros(ndd,size(ss,2)); % Transmitted symbol buffer
MaxIter=4; % # of iterations for getting the probability of error
nobe= 0; % Number of bit errors to be accumulated
thh(1)=0; nh(1)=-1; % Initial estimate of phase offset/sampling instant
ycsk0=[0  0]; % Previous correlator output samples initialized
rand('twister',5489); randn('state',0); 
for k=1:MaxIter
i= ceil(rand*M); s=ss(i,:); % Signal index, Data bits to transmit
[ch_input,zi_XMTR]= filter(g,1,sw(i,:),zi_XMTR);
wct=wc*T*[0:Ns-1];
noise= sgmsT*(randn(1,Ns).*cos(wct)-randn(1,Ns).*sin(wct));
r = ch_input + noise; % Received signal
yr=[yr(:,Ns+1:LB) suT.*[r;r]]; % Multiplier output (buffer)
ycsk=sum(yr(:,interval+nh(k))'); % Correlator output at (k-1)Ts
ycsk1=sum(yr(:,interval1+nh(k))'); % Correlator output at (k-1/2)Ts
e(k)=ycsk*(ycsk1-ycsk0).';  % Eq.(8.5.3a)
dnh=0; % Stepsize of adjusting the sampling instant
if KC>1 % For timing recovery
     if e(k)>eth, dnh=1; elseif e(k)<-eth, dnh=-1; end 
end  
nh(k+1)=nh(k)+dnh; % Eq.(8.5.2) for timing recovery
ycsk0=ycsk1; % Correlator output sample buffer updated
ynejth= ycsk(1)+j*ycsk(2);
if KC>0, ynejth=ynejth*exp(-j*thh(k)); end % Phase compensation
th=angle(ynejth); % Compensated phase of received signal waveform
if th<-pi/M, th=th+2*pi; end
[themin,imin]=min(abs(th-phases)); D=ss(imin,:); % Detected data 
thh(k+1)=thh(k)+gamma*imag(exp(-j*phases(imin))*ynejth); %(8.4.10)
if k>ndd,  nobe=nobe+sum(s0s(1,:)~=D);  end 
if nobe>100, break; end
s0s=[s0s(2:end,:); s]; % Transmitted symbol buffer updated
end
pobe= nobe/((k-ndd)*b) % The BER