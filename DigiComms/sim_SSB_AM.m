function sim_SSB_AM(U_or_L)
% simulates USSB/LSSB-AM depending on the value of USSB_or_LSSB=0/1
if nargin<1, U_or_L=0;  end
clear, clf
Ac=1; fc=50; wc=2*pi*fc;     % Amplitude/Frequency of carrier
Tb=0.1;                      % Bit interval time
T=1/fc/8; Fs=1/T;            % Sampling period/frequency
Nb=Tb/T; lt=2^(nextpow2(3*Nb)); t=[1:lt]*T; % Time vector
m= ones(Nb,1)*[4 -8 -4]; m=m(:).';          % Message signal m(t)
m=[m, zeros(1,lt-length(m))];  ma=hilbert(m);  
tmpc=real(ma).*cos(wc*t); tmps=imag(ma).*sin(wc*t);
if U_or_L<=0;  m_ssb=Ac*(tmpc-tmps); str='USSB-AM'; % USSB-AM signal
 else     m_ssb=Ac*(tmpc+tmps); str='LSSB-AM';  % LSSB-AM signal
end
y_ssb=m_ssb*2/Ac.*cos(wc*t); % Demodulated signal
% Digital FIR LPF design
Bd= fir1(20,fc*T); Ad=1;
% Output of LPF as a detector
y_dtr=filter(Bd,Ad,y_ssb); 
clf, plot_MOD(T,lt,m,m_ssb,y_ssb,str,y_dtr)

%t=[1:lt]*T;
m_ssb=modulate(m,fc,Fs,'amssb');
subplot(423), hold on, plot(t,m_ssb,'r')
y_ssb=demod(m_ssb,fc,Fs,'amssb');
subplot(427), hold on, plot(t,y_ssb,'r')



