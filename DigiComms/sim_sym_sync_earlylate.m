%sim_sym_sync_earlylate.m -- Early-Late sampling time control
clear, clf
KC=0; % 1/0 for timing recovery or not
Ns=80; % Number of sample times per symbol time
Ts=1; T=Ts/Ns; % Symbol/Sample time
Delay=3; % Number of delay symbol times of XMTR filter
rof=0.5;  % Roll-off factor
% Square Root Raised Cosine Filtering with upsampling of Ns ? P6.1(h)
Fd=1/Ts; Fs=Fd*Ns; SRRC_filter = rcosine(Fd,Fs,'fir/sqrt',rof,Delay);
% SRRC Filtering with upsampling 
[g,tg]=rcosflt([1],Fd,Fs,'filter',SRRC_filter); 
g=g.'; Ng=length(g); % Channel response and its length
h=fliplr(g); % Impulse response of RCVR matched filter in baseband
Es=1; b=1; SNRbdB=8; b=1; pobet= prob_error(SNRbdB,'PSK',b,'bit')
SNRb=10^(SNRbdB/10); N0=2*(Es/b)/SNRb; sigma2=N0/2; sgmsT=sqrt(sigma2/T);
dT=Ts/2;  dTN=ceil(dT/T);  delta=0.01;
MaxIter=2000;
for ts=[3.1 2.9 2.1]*Ts
   tsN=ceil(ts/T); % Initial sampling time in number of sample times
   dd=round(Ng/Ns)*2-floor(ts/Ts)+2; % Detection delay
   lch=Ng+Ns-1; ch_input=zeros(1,lch); % Channel input buffer
   lrs=Ng+Ns-1; rs=zeros(1,lrs); lys=Ns*10; ys=zeros(1,lys);
   lds=100+dd;  ds=zeros(1,lds); lDs=100;   Ds=zeros(1,lDs);
   nobe=0; 
   for k=1:MaxIter
      d(k)=(rand>0.5)*2-1; % BPSK-modulated data symbol
      ds=[ds(2:end) d(k)];
      ch_input= [ch_input(Ns+1:end) d(k)*ones(1,Ns)];
      for n=1:Ns % through channel
         rn=ch_input(end-Ns+n-Ng+1:end-Ns+n)*g' +sgmsT*randn;
         rs=[rs(2:end) rn/Ns];  % Received signal buffer updated
         yn=rs(end-Ng+1:end)*h'; % RCVR Matched filter output
         ys=[ys(2:end) yn]; % RCVR filter output buffer updated
      end
      dif=abs(ys(tsN-dTN))-abs(ys(tsN+dTN));
      if KC>0 % Early or Late by Eq.(8.5.1)
        if dif>delta, ts=ts-T; elseif dif<-delta, ts=ts+T; end
      end  
      ts_history(k)=ts; tsN=ceil(ts/T); % Adjusted sampling time
      D=(ys(tsN)>0)*2-1; Ds= [Ds(2:lDs) D]; % Detection and DTR buffer nobe=nobe+(D~=ds(lds-dd)); if nobe>100, break; end
   end  
   pobe=nobe/k
   plot(ts_history) 
end 