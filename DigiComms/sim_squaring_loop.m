%sim_squaring_loop.m
% simulates a PLL with squaring loop (Fig. 8.5)
clear, clf
Ts=1e-5; Ns=64; T=Ts/Ns;
wc=8*pi/Ts;  
ds=randint(1,20)*2-1; % A sequence of BPSK symbols
N_ds=length(ds); NT= N_ds*Ns;
th=pi/5 +pi/20*sin([1:NT]*pi/200); % Time-varying phase of input 
% Elliptic Digital Filter Design
fp=0.05; fs=0.1; Rp=3; As=40;
[Norder,fc0]= ellipord(fp,fs,Rp,As); % Order of digital elliptic filter
[Bde,Ade]= ellip(Norder,Rp,As,fc0); % Num/Den of elliptic filter
NB=length(Bde); NA=length(Ade);
lb_ve=10; b_ve=zeros(1,lb_ve); % initialize buffer for ve
thh_s2(1)=0;
vLF(1)=0; KvT=10/Ts*T; 
for KC=0:1
   n=1; 
   zi=zeros(1,max(NA,NB)-1); % Initial state of the LPF
   for k=1:N_ds  % Symbol-by-symbol processing
      d(k)=ds(k); % BPSK symbol to transmit
      for m=1:Ns
         t(n)=(n-1)*T;  wct=wc*t(n); wct2=wct*2;
         dt(n)=d(k); r(n)= d(k)*cos(wct+th(n)); % Received signal
         if KC==0 % multiplying with non-phase-compensated carrier
           r_dem(n)=r(n)*2*cos(wct); 
          else % multiplying with phase-compensated carrier
           vVCO(n)= cos(wct2+thh_s2(n)); % VCO output
           ve(n)= r(n)^2*vVCO(n); % Input to LF (Eq.(8.4.2))
           vLF(n+1)= vLF(n)+(ve(n)-b_ve(1))*KvT; %LF output
           b_ve= [b_ve(2:lb_ve) ve(n)]; % Buffer for ve(t)
           thh_s2(n+1)= thh_s2(n)+KvT*vLF(n+1); % Phase estimate
           % r(t) multiplied by the frequency divider output
           thh(n)= (thh_s2(n)-pi/2)/2; % Eq.(8.4.5)
           r_dem(n)=r(n)*2*cos(wct+thh(n)); % Demodulation
         end 
         [dh(n),zi]=filter(Bde,Ade,r_dem(n),zi); % sample by sample
         n=n+1;
      end     
   end
   if KC>0,  subplot(311), plot(t,th,'k', t,thh,'r');  end
   subplot(312+KC), plot(t,dt,'k', t,dh,'r')
end