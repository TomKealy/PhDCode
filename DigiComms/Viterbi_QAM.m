function [pemb,nombe,notmb]=Viterbi_QAM(Gc,b,SNRbdB,MaxIter)
if nargin<4, MaxIter=1e5; end
if nargin<3, SNRbdB=5; end
if nargin<2, b=4; end
[K,N]=size(Gc); Rc=K/N;  Gc_m=max(Gc.');
% Constraint length vector 
for i=1:length(Gc_m), Lc(i)=length(deci2bin1(oct2dec(Gc_m(i))));
Nf=144; % Number of bits per frame
Nmod=Nf*N/K/b; % Number of QAM symbols per modulated frame
SNRb=10.^(SNRbdB/10); SNRbc=SNRb*Rc; sqrtSNRbc=sqrt(SNRbc); 
sqrtSNRc=sqrt(2*b*SNRbc); % Complex noise for b-bit (coded) symbol
trel=poly2trellis(Lc,Gc); 
Tbdepth=5; delay=Tbdepth*K;
nombe=0; Target_no_of_error=100;
for iter=1:MaxIter
   msg=randint(1,Nf); % Message vector
   coded= convenc(msg,trel); % Convolutional encoding
   modulated= QAM(coded,b); % 2^b-QAM-Modulation
   r= modulated +(randn(1,Nmod)+j*randn(1,Nmod))/sqrtSNRc;
   demodulated= QAM_dem(r,b); % 2^b-QAM-Demodulation
   decoded= vitdec(demodulated,trel,Tbdepth,'trunc','hard');
   nombe = nombe + sum(msg~=decoded(1:Nf));  %
if nombe>Target_no_of_error, break; end
end
notmb=Nf*iter; % Number of total message bits
pemb=nombe/notmb; % Message bit error probability