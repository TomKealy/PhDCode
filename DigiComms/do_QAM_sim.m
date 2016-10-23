%do_QAM_sim.m
clear, clf
b=4; M=2^b; L=2^(b/2); % # of bits per symbol and Modulation order
Ns=b*10; Ts=1e-5; T=Ts/Ns; % Symbol time and Sample time
wc=2*pi*10/Ts; % Carrier Frequency[rad/s] (such that wc*T<pi)
Target_no_of_error=100; 
SNRdBs=[5  10]; EbN0dBs=SNRdBs-3;
for i=1:length(SNRdBs)
SNRdB=SNRdBs(i); EbN0dB=SNRdB-3; % Eb/(N0/2)=SNR-> Eb/N0=SNR/2
sim('QAM_passband_sim',1e5*Ts); SERs(i) = ser(end,1);
end
SNRdBt=0:0.1:13; poset = prob_error(SNRdBt,'QAM',b,'sym'); 
semilogy(SNRdBt,poset,'k', SNRdBs,SERs,'*'), xlabel('SNR[dB]')
Transmitted_Signal_Power = 2*(M-1)/3/Ts;
Received_Signal_Power = Transmitted_Signal_Power*(1+10^(-(EbN0dBs(end)+3)/10)*Ts/b/T)
