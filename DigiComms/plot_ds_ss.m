%plot_ds_ss.m (to draw Fig. 10.5)
clear, clf
c=[1 0 1]; % A simple PN code, i.e., c=png([1,4]); 
Nc=length(c);
T=0.1/16; Tc=0.1; NTc=Tc/T; Tb=Nc*Tc; NTb=Nc*NTc; % Chip/Bit interval
%d1=ones(1,NTc); d1b=ones(1,NTb); % For oversampling
d=[1 0 1 1]; d_bp=2*d-1; % Channel coded message & its bipolar version
dt=ones(Nc*NTc,1)*d; dt=dt(:)'; % Data signal d(t)
d_bpt=2*dt-1; % Bipolar/Non-Return_to_Zero data signal
N=length(dt); % 4*NTb=4Tb/T
t=[0:4*NTb-1]*T; tf=4*NTb*T;
subplot(621), plot(t,d_bpt,'k'), axis([0 tf -1.3 1.3])
f=[-N/2:N/2]*2*pi/N/T;
Df= fftshift(abs(fft(d_bpt,N))); Df=[Df Df(1)];
subplot(622), plot(f,Df,'k')
c_bp=2*c-1; c_bpt=ones(NTc,1)*c_bp; % Bipolarized code
dct=c_bpt(:)*d_bp; dct=dct(:)'; % d(t)c(t): PN-coded seq.
subplot(623), plot(t,dct,'k'), axis([0 tf -1.3 1.3])
Dcf= fftshift(abs(fft(dct,N))); Dcf=[Dcf Dcf(1)];
subplot(624), plot(f,Dcf,'k')
wc=60*pi; % Carrier frequency
carrier=cos(wc*[0:NTc-1]*T);
dcmt=carrier'*dct([NTc:NTc:4*NTb]); dcmt=dcmt(:)'; % Modulated
subplot(625), plot(t,dcmt,'k'), axis([0 tf -1.3 1.3])
Dcmf= fftshift(abs(fft(dcmt,N))); Dcmf=[Dcmf Dcmf(1)];
subplot(626), plot(f,Dcmf,'k')
dcmdt=xcorr1(dcmt,2*carrier); % Demodulated
dcmdd= dcmdt([NTc:NTc:4*NTb]); % Sampled
dcmddt=ones(NTc,1)*dcmdd; dcmddt=dcmddt(:)';
subplot(627), plot(t,dcmddt,'k'), axis([0 tf -1.3 1.3])
Dcmddf= fftshift(abs(fft(dcmddt,N))); Dcmddf=[Dcmddf Dcmddf(1)];
subplot(628), plot(f,Dcmddf,'k')
dcc=c_bp*reshape(dcmdd,Nc,length(dcmdd)/Nc); % Decoded
dccmdt=ones(NTb,1)*dcc(:)'/Nc; 
dccmdt=dccmdt(:)'; % Oversampled and decoded
subplot(629), plot(t,dccmdt,'k'), axis([0 tf -1.3 1.3])
Dccmdf= fftshift(abs(fft(dccmdt,N))); Dccmdf=[Dccmdf Dccmdf(1)];
subplot(6,2,10), plot(f,Dccmdf,'k')
