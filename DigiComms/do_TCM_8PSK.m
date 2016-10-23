%do_TCM_8PSK.m
clear 
Nf=144; Tf=0.001; Tb=Tf/Nf; % Frame size, Frame and Sample/Bit times
Gc=[1 0 0; 0 5 2]; state_eq='TCM_state_eq1'; Ns=2; % Fig. P9.9.1(a)
%Gc=[1 2 0; 4 1 2]; state_eq='TCM_state_eq2'; Ns=3; % Fig. P9.9.1(b)
[K,N]=size(Gc); Rc=K/N; % Message/Codeword length and Code rate
% Constellation for TCM Simuklink block
Constellation2=exp(j*2*pi/8*[0 4 2 6 1 5 3 7]); 
% Constellation good for TCM() MATLAB function
Constellation1=exp(j*2*pi/8*[0:7]); 
% Constraint length vector 
Gc_m=max(Gc.');
for i=1:length(Gc_m), Lc(i)=length(deci2bin(oct2dec(Gc_m(i))));  end
%trel=poly2trellis(Lc,Gc);
Tbdepth=sum(Lc)*5; delay=Tbdepth*K;
Ts=K*Tb; % or N*Rc*Tb % Symbol time
EbN0dBs_t=0:0.1:10;  SNRbdBs_t=EbN0dBs_t+3;
BER_theory= prob_error(SNRbdBs_t,'PSK',K,'BER');
EbN0dBs=2; %[3 6];
Target_no_of_error=50;
for i=1:length(EbN0dBs)
   EbN0dB=EbN0dBs(i); SNRbdB=EbN0dB+3; 
   pemb=TCM(state_eq,K,Ns,N,Constellation1,SNRbdB,Target_no_of_error);
   pembs_TCM8PSK(i)=pemb
   sim('TCM_sim'); pembs_TCM8PSK_sim(i)=BER(1);
end
figure(1), clf
semilogy(EbN0dBs_t,BER_theory,'b'), hold on
semilogy(EbN0dBs,pembs_TCM8PSK_sim,'r*', EbN0dBs,pembs_TCM8PSK,'m*')
