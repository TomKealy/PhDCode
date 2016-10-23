%sim_TCM.m 
% simulates a Trellis-coded Modulation
clear, clf
lm=1e4; msg=randint(1,lm);
K=2; N=K+1; Ns=3; % Size of encoder input/output/state
M=2^N; Constellation=exp(j*2*pi/M*[0:M-1]); % Constellation
ch_input=TCM_encoder('TCM_state_eq0',K,Ns,N,msg,Constellation);
lc= length(ch_input);
SNRbdB=5; SNRdB=SNRbdB+10*log10(K); % SNR per K-bit symbol
sigma=1/sqrt(10^(SNRdB/10)); noise=sigma/sqrt(2)*(randn(1,lc)+j*randn(1,lc)); % Complex noise
received = ch_input + noise; var(received)
received = awgn(ch_input,SNRdB); var(received) % Alternative
decoded_seq=TCM_decoder('TCM_state_eq0',K,Ns,received,Constellation);
ber_TCM8PSK = sum(decoded_seq(1:lm)~=msg)/lm
ber_QPSK_theory = prob_error(SNRbdB,'PSK',K,'BER')