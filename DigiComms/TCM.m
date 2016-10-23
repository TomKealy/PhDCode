function [pemb,nombe,notmb]=TCM(state_eq,K,Nsb,N,Constellation,SNRbdB,Target_no_of_error)
% Copyleft: Won Y. Yang, wyyang53@hanmail.net, CAU for academic use only
M=2^N; Rc=K/N;
Nf=288; % Number of bits per frame
Nmod=Nf/K; % Number of symbols per modulated frame
SNRb=10.^(SNRbdB/10); SNRbc=SNRb*Rc; 
sqrtSNRc=sqrt(2*N*SNRbc); % Complex noise per K=Rc*N-bit symbol
nombe=0;  MaxIter=1e6; 
for iter=1:MaxIter
   msg=randint(1,Nf); % Message vector
   coded = TCM_encoder(state_eq,K,Nsb,N,msg,Constellation); 
   r= coded +(randn(1,Nmod)+j*randn(1,Nmod))/sqrtSNRc; 
   decoded= TCM_decoder(state_eq,K,Nsb,r,Constellation);
   nombe = nombe + sum(msg~=decoded(1:Nf));
   if nombe>Target_no_of_error, break; end
end
notmb=Nf*iter; % Number of total message bits
pemb=nombe/notmb; % Message bit error probability