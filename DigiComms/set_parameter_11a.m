%set_parameter_11a.m
% sets the parameters of IEEE 802.11a with Synchronization
clear all;
EbN0dB = 10;  % Eb/N0(dB)
CFO = 1.2;  % An imagined CFO
STO = 0;  % Symbol timing offset
Nfft=64; Ng=16; Nsym=Ng+Nfft;  % FFT size, Guard size, OFDM symbol size
Tsym=4e-6; T=Tsym/Nsym; % OFDM Symbol Time and Sample Time
% Null symbol, short preamble, and long preamble
null_symbol = zeros(1,Nsym)+eps*i;
[short_preamble,Short]=short_train_seq(Nfft);
[long_preamble,Long]=long_train_seq(Nfft);
Data_rate=36,  Nbpsc=1; Nsd=48; Ncbps = Nsd*Nbpsc; 
load ch_80.dat
ch=ch_80([1:12],1)+ch_80([1:12],2)*i;
H_true = fft(ch,Nfft);  H_true=H_true([39:64 2:27]);