%start_PSK_carrier_phase_timing_recovery.m
clear
b=2; M=2^b; SNRbdB=10;
%b=3; M=2^b; SNRbdB=12;
Ts=1e-5; Ns=8; Nf=100;
time_delay=1.5; % in number of sample
phase_offset=10; % in degree
SNRdB=SNRbdB+10*log10(b/Ns)-3
EsN0dB=SNRbdB+10*log10(b)-3
EbN0dB=SNRbdB-3
pobet=prob_error(SNRbdB,'PSK',b,'sym')
rof=0.4; SRC_Delay=4;
sim('PSK_carrier_phase_timing_recovery',1e6*Ts)
