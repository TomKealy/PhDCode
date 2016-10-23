%do_FSK_sim.m
clear, clf
b=1; M=2^b; % Number of bits per symbol and Modulation order
Nbsps=b*2^3;  % # of samples per symbol in baseband
Nos=5; Ns=Nbsps*Nos; % # of samples per symbol in passband
Ts=1e-5; T=Ts/Ns; % Symbol time and Sample time
dw=2*pi/Ts; % Frequency separation [rad/s]
wc=10*dw; % Carrier Frequency[rad/s] (such that wc*T<pi)
Target_no_of_error=50; SNRdBs=[5  10]; EbN0dBs=SNRdBs-3;
for i=1:length(SNRdBs)
   SNRdB=SNRdBs(i); EbN0dB=SNRdB-3; % Eb/(N0/2)=SNR-> Eb/N0=SNR/2
   sim('FSK_passband_sim',1e5*Ts); SERs(i)=ser(end,1);
end
SNRdBt=0:0.1:13; 
poset = prob_error(SNRdBt,'FSK',b,'sym','noncoherent'); 
semilogy(SNRdBt,poset,'k', SNRdBs,SERs,'*'), xlabel('SNR[dB]')
Transmitted_Signal_Power = 1/Ts;
Received_Signal_Power = Transmitted_Signal_Power*(1+10^(-(EbN0dBs(end)+3)/10)*Ts/b/T)
