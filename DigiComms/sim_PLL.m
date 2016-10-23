%sim_PLL.m 
% simulates a PLL with loop filter
clear, clf
Ts=1e-5; Ns=40; T=Ts/Ns; % Symbol time and Sample time
MaxIter=400; tt=[0:MaxIter]*T; % Time range of simulation
wc=8*pi/Ts; pi2=pi*2;
% Changing phase of input 
th= pi/4 + pi/50*sin(0.01*pi/T*tt); 
% Desired damping ratio & natural freq of PLL
zeta=0.707; wn=pi/2/Ts; 
% Parameters of Loop Filter to be designed
T2=max(2*zeta/wn,T/pi), T1=10*T2 % Eq.(8.2.9,10)
K=2*wn^2*T1;  % Eq.(8.2.9)
% Discretize the LF (T2*s+1)/(T1*s) through BLT with prewarping
[bLF,aLF]=bilinear([T2 1],[T1 0],1/T,1/T2/2/pi)
% Discretize the VCO (Kv/s) through BLT without prewarping
[bVCO,aVCO]=bilinear(K,[1 0],1/T)
% Initialize the phase estimate & integral of LF output
thh(1)=0; vLF(1)=0; 
for n=1:MaxIter
   t=tt(n);  vi(n)=sin(wc*t+th(n)); % Input to PLL 
   vVCO(n)=cos(wc*t+thh(n)); % VCO output 
   ve(n)=vi(n)*vVCO(n);     % Input to LF (Eq.(8.2.4))
   if n<=1
     vLF(n+1)=-aLF(2)*vLF(n)+bLF(1)*ve(n); % LF output (Eq.(8.2.10))
   else  
     vLF(n+1)=-aLF(2)*vLF(n)+bLF*ve(n:-1:n-length(bLF)+1)'; 
   end    
   % Phase estimate made by VCO 
   thh(n+1)= -aVCO(2)*thh(n)+bVCO*vLF(n+1:-1:n-length(bVCO)+2)'; 
   if thh(n+1)>pi,  thh(n+1)= thh(n+1)-pi2;
    elseif thh(n+1)<-pi,  thh(n+1)= thh(n+1)+pi2;
   end
end
subplot(211), plot(tt,thh,'r', tt,th,'k')

% To simulate a PLL with loop filter replaced by an integrator
lb_ve=10; b_ve=zeros(1,lb_ve); % Initialize the buffer for ve(t)
thh(1)=0; vLF(1)=0; KT=1/Ts;
for n=1:MaxIter
   t=tt(n);  vi(n)=sin(wc*t+th(n)); % Input to PLL
   vVCO(n)=cos(wc*t+thh(n)); % VCO output
   ve(n)= vi(n)*vVCO(n); % Input to LF
   vLF(n+1)=vLF(n)+(ve(n)-b_ve(1))*T; % LF output (Eq.(8.2.12))
   b_ve=[b_ve(2:lb_ve) ve(n)]; % Buffer for ve
   thh(n+1)= thh(n)+KT*vLF(n+1); % Phase estimate (Eq.(8.2.11))
   if thh(n+1)>pi,  thh(n+1)= thh(n+1)-2*pi;
    elseif thh(n+1)<-pi,  thh(n+1)= thh(n+1)+2*pi;
   end
end
subplot(212), plot(tt,th,'k', tt,thh,'r')