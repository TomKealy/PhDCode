%do_BCH_BPSK_sim.m
clear, clf
K=16; % Number of input bits to the BCH encoder (message length)
N=31; % Number of output bits from the BCH encoder (codeword length)
Rc=K/N; % Code rate to be multiplied with the SNR in AWGN channel block
b=1; M=2^b; % Number of bits per symbol and modulation order
T=0.001/K; Ts=b*T; % Sample time and Symbol time
SNRbdBs=[2:4:10]; EbN0dBs=SNRbdBs-3;
EbN0dBs_t=0:0.1:10; EbN0s_t=10.^(EbN0dBs_t/10); 
SNRbdBs_t=EbN0dBs_t+3;
BER_theory= prob_error(SNRbdBs_t,'PSK',b,'BER');
for i=1:length(EbN0dBs)
   EbN0dB=EbN0dBs(i);
   sim('BCH_BPSK_sim');
   BERs(i)=BER(1); % just ber among {ber, # of errors, total # of bits}
   fprintf(' With EbN0dB=%4.1f, BER=%10.4e=%d/%d\n', EbN0dB,BER);
end
semilogy(EbN0dBs,BERs,'r*', EbN0dBs_t,BER_theory,'b')
xlabel('Eb/N0[dB]'); ylabel('BER'); 
title('BER of BCH code with BPSK');