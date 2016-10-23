%do_Viterbi_QAM.m
clear, clf
Nf=144; Tf=0.001; Tb=Tf/Nf; % Frame size, Frame and Sample/Bit time
Gc=[133 171]; 
[K,N]=size(Gc); Rc=K/N; % Message/Codeword length and Code rate
% Constraint length vector 
Gc_m=max(Gc.');
for i=1:length(Gc_m)
Lc(i)=length(deci2bin1(oct2dec(Gc_m(i)))); 
end
Tbdepth=sum(Lc)*5; delay=Tbdepth*K;
b=4; M=2^b; % Number of bits per symbol and Modulation order
Ts=b*Rc*Tb; % Symbol time
N_factor=sqrt(2*(M-1)/3); % Eq.(7.5.4a)
EbN0dBs_t=0:0.1:10; SNRbdBs_t=EbN0dBs_t+3;
BER_theory= prob_error(SNRbdBs_t,'QAM',b,'BER');
EbN0dBs=[3 6]; Target_no_of_error=50;
for i=1:length(EbN0dBs)
   EbN0dB=EbN0dBs(i); SNRbdB=EbN0dB+3; 
   randn('state', 0);
   [pemb,nombe,notmb]=???????_QAM(Gc,b,SNRbdB,Target_no_of_error); pembs(i)=pemb;
   sim('Viterbi_QAM_sim'); pembs_sim(i)=BER(1); 
end
[pembs; pembs_sim] 
semilogy(EbN0dBs,pembs,'r*', EbN0dBs_t,BER_theory,'b')
xlabel('Eb/N0[dB]'); ylabel('BER');